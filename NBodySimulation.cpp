#include "NBodySimulation.h"

NBodySimulation::NBodySimulation () :
  t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
  x(nullptr), v(nullptr), mass(nullptr),
  timeStepSize(0), maxV(0), minDx(0), videoFile(nullptr),
  
  snapshotCounter(0), timeStepCounter(0) {};

NBodySimulation::~NBodySimulation () {
  if (x != nullptr) {
    for (int i=0; i<NumberOfBodies; i++)
      delete [] x[i];
    delete [] x;
  }
  if (v != nullptr) {
    for (int i=0; i<NumberOfBodies; i++)
      delete [] v[i];
    delete [] v;
  }
  if (mass != nullptr) {
    delete [] mass;
  }

}

void NBodySimulation::checkInput(int argc, char** argv) {
    if (argc==1) {
    std::cerr << "usage: " << std::string(argv[0])
              << " plot-time final-time dt objects" << std::endl
              << " Details:" << std::endl
              << " ----------------------------------" << std::endl
              << "  plot-time:       interval after how many time units to plot."
                 " Use 0 to switch off plotting" << std::endl
              << "  final-time:      simulated time (greater 0)" << std::endl
              << "  dt:              time step size (greater 0)" << std::endl
              << "  objects:         any number of bodies, specified by position, velocity, mass" << std::endl
              << std::endl
              << "Examples of arguments:" << std::endl
              << "+ One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "    0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0" << std::endl
              << "+ One body spiralling around the other" << std::endl
              << "    0.05  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0" << std::endl
              << "+ Three-body setup from first lecture" << std::endl
              << "    0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0" << std::endl
              << "+ Five-body setup" << std::endl
              << "    0.02  20.0  0.0009    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0" << std::endl
              << std::endl;

    throw -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each body is given by seven entries"
                 " (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments"
                 " (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    throw -2;
  }
}

void NBodySimulation::setUp (int argc, char** argv) {

  checkInput(argc, argv);

  NumberOfBodies = (argc-4) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies"
            << std::endl;

  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta
              << " time units" << std::endl;
    tPlot = 0.0;
  }
}

double NBodySimulation::euclidean_distance(int i,int j){
    const double distance = sqrt(
                               (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
                               (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
                               (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
                               );
    return distance;
}

// on_particle_idx: index of the particle the force is acting ON.
// by_particle_idx: index of the particle EXERTING the force.
// direction: 0 for x, 1 for y, 2 for z.
// precomputed_distance: the pre-calculated distance between on_particle_idx and by_particle_idx.
double NBodySimulation::force_calculation(int on_particle_idx, int by_particle_idx, int direction, double precomputed_distance) {
    if (precomputed_distance == 0.0) {
        return 0.0; // Avoid division by zero if distance is zero.
    }
    double distance3 = precomputed_distance * precomputed_distance * precomputed_distance;
    // Handle cases where distance is extremely small, leading to distance3 being zero.
    if (distance3 == 0.0) { 
         return 0.0;
    }

    // The force ON on_particle_idx BY by_particle_idx has components proportional to 
    // (x[by_particle_idx][direction] - x[on_particle_idx][direction]).
    // The original code's force_calculation(j,i,dir) returned (x[j][dir]-x[i][dir])*m[j]*m[i]/d^3.
    // To maintain this, if called as force_calculation(j, i, dir, dist),
    // on_particle_idx = j, by_particle_idx = i.
    // We need to return (x[on_particle_idx][direction] - x[by_particle_idx][direction]) * mass[on_particle_idx] * mass[by_particle_idx] / distance3
    // to keep the existing summation logic in updateBody correct.
    // (x[j_orig][dir] - x[i_orig][dir]) * mass[j_orig] * mass[i_orig] / d^3
    return (x[on_particle_idx][direction] - x[by_particle_idx][direction]) * mass[on_particle_idx] * mass[by_particle_idx] / distance3;
}

double NBodySimulation::NBodySimulationMolecularForces(int i, int j, int direction) {
    double rij[3]; 
    double rij_norm_sq = 0.0;

    for (int d = 0; d < 3; d++) {
        rij[d] = x[i][d] - x[j][d];
        rij_norm_sq += rij[d] * rij[d];
    }

    if (rij_norm_sq == 0.0) { // If distance is zero
        return 0.0;
    }
    double rij_norm = sqrt(rij_norm_sq);

    // Denominators for the force formula terms
    double d_pow9 = pow(rij_norm, 9);
    double d_pow13 = pow(rij_norm, 13);

    // If rij_norm is so small that powers result in zero (underflow) or infinity (overflow, less likely for positive powers)
    // or if rij_norm itself is zero (already handled by rij_norm_sq check for exact zero),
    // prevent division by zero.
    if (d_pow9 == 0.0 || d_pow13 == 0.0) {
        // This implies rij_norm is extremely small. The force would be theoretically huge.
        // Returning 0.0 to avoid NaN/inf, though a large capped value might be more "physical"
        // depending on the simulation's needs for such close encounters.
        return 0.0; 
    }

    double force_magnitude = 10.0 * ((0.1 / d_pow13) - (0.1 / d_pow9)) * rij_norm;
    
    // The final division by rij_norm for component is safe because rij_norm is not zero here.
    double force_component = force_magnitude * rij[direction] / rij_norm;

    return force_component;
}

void NBodySimulation::updateBody () {
  // std::cout << "energy: "<< calculate_energy() <<std::endl;
  timeStepCounter++;
  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();

  double half_dt = timeStepSize/2;

// compute half an Euler time step for dV --------------------------------------
  double* force0 = new double[NumberOfBodies]();
  double* force1 = new double[NumberOfBodies]();
  double* force2 = new double[NumberOfBodies]();

    for (int i = 0; i < NumberOfBodies; i++) {
      for (int j = i+1; j < NumberOfBodies; j++) {
          
          const double distance = this->euclidean_distance(i, j);
          minDx = std::min(minDx, distance);

          double Fx = force_calculation(j, i, 0, distance);
          double Fy = force_calculation(j, i, 1, distance);
          double Fz = force_calculation(j, i, 2, distance);


          force0[i] += Fx;
          force1[i] += Fy;
          force2[i] += Fz;


          force0[j] -= Fx;
          force1[j] -= Fy;
          force2[j] -= Fz;
      }
  }

  for (int i = 0; i<NumberOfBodies;i++){
    v[i][0] = v[i][0] + half_dt * force0[i] / mass[i];
    v[i][1] = v[i][1] + half_dt * force1[i] / mass[i];
    v[i][2] = v[i][2] + half_dt * force2[i] / mass[i];
  }

//update positions  --------------------------------------
  for (int i = 0; i<NumberOfBodies;i++){
    x[i][0] = x[i][0] + timeStepSize * v[i][0];
    x[i][1] = x[i][1] + timeStepSize * v[i][1];
    x[i][2] = x[i][2] + timeStepSize * v[i][2];
  }
// update force from new position --------------------------------------
  for (int i = 0; i < NumberOfBodies; i++) {
      force0[i] = 0.0;
      force1[i] = 0.0;
      force2[i] = 0.0;
  }

  for (int i = 0; i < NumberOfBodies; i++) {
      for (int j = i+1; j < NumberOfBodies; j++) {

          const double distance = this->euclidean_distance(i, j);
          minDx = std::min(minDx, distance);

          double Fx = force_calculation(j, i, 0, distance);
          double Fy = force_calculation(j, i, 1, distance);
          double Fz = force_calculation(j, i, 2, distance);


          force0[i] += Fx;
          force1[i] += Fy;
          force2[i] += Fz;


          force0[j] -= Fx;
          force1[j] -= Fy;
          force2[j] -= Fz;
      }
  }

//update the velocities --------------------------------------
  for (int i = 0; i<NumberOfBodies;i++){
    v[i][0] = v[i][0] + half_dt * force0[i] / mass[i];
    v[i][1] = v[i][1] + half_dt * force1[i] / mass[i];
    v[i][2] = v[i][2] + half_dt * force2[i] / mass[i];
    double cur_v = std::sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2] );
    
    maxV = std::max(cur_v,maxV);

  }

  //coll
   // C_factor calculation as before, ensuring it's defined based on current NumberOfBodies.
   double C_factor = (NumberOfBodies > 0) ? 0.01/NumberOfBodies : 0.0;

    for (int i = 0; i < NumberOfBodies; i++) {
      for (int j = i+1; j < NumberOfBodies; j++) {
          // Safety break if all bodies somehow get eliminated (though unlikely with this logic)
          if (NumberOfBodies == 0) break; 

          double threshold_distance = C_factor * (mass[i] + mass[j]);

          if (euclidean_distance(i,j) <= threshold_distance){
            std::cout << "Collision between body " << i << " and body " << j << std::endl;
            
            // Calculate new properties for body i (merged body)
            // Important: Use the state of mass[i] *before* it's updated for the weighted average.
            double old_mass_i = mass[i];
            double combined_mass = old_mass_i + mass[j];
            
            if (combined_mass == 0) { 
                std::cerr << "Error: Combined mass is zero during collision between " << i << " and " << j << ". Skipping collision." << std::endl;
                continue; 
            }

            for (int direction = 0; direction <= 2; direction++) {
                x[i][direction] = (old_mass_i*x[i][direction] + mass[j]*x[j][direction]) / combined_mass;
                v[i][direction] = (old_mass_i*v[i][direction] + mass[j]*v[j][direction]) / combined_mass;
            }
            mass[i] = combined_mass; // Update mass of body i

            // Now handle removal of body j
            // 1. Free memory of the particle that is disappearing (original body j)
            delete[] x[j];
            delete[] v[j];

            int lastBodyIndex = NumberOfBodies - 1; // This is the index of the last particle *before* decrementing NumberOfBodies

            // 2. If j is not the last particle, move data from last particle to slot j
            if (j < lastBodyIndex) {
                x[j] = x[lastBodyIndex];
                v[j] = v[lastBodyIndex];
                mass[j] = mass[lastBodyIndex];
                // The pointers x[lastBodyIndex] and v[lastBodyIndex] are now effectively stale
                // as their content has been moved. They will be outside the valid range
                // after NumberOfBodies is decremented. No need to null them.
            }
            // If j == lastBodyIndex, its memory was freed, and no data needs to be moved.
            // The slot x[lastBodyIndex], v[lastBodyIndex] is now gone.

            // 3. Decrement number of bodies
            NumberOfBodies--;

            // 4. Adjust loop variable j
            // This ensures that if a body was moved into slot j, it gets checked against body i.
            // It also correctly adjusts for the reduced size of NumberOfBodies for subsequent iterations.
            j--; 
          }
      }
  }
  // std::cout << std::sqrt( v[0][0]*v[0][0] + v[0][1]*v[0][1] + v[0][2]*v[0][2] ) <<std::endl;
  // if (timeStepCounter  % 100000 == 0 ){
  //   std::cout << "energy: "<< calculate_energy() <<std::endl;
  //   calculate_momentum();
  //   }
  t += timeStepSize;
  delete[] force0;
  delete[] force1;
  delete[] force2;
}

/**
 * Check if simulation has been completed.
 */
bool NBodySimulation::hasReachedEnd () {
  return t > tFinal;
}

void NBodySimulation::takeSnapshot () {
  if (t >= tPlot) {
    printParaviewSnapshot();
    printSnapshotSummary();
    tPlot += tPlotDelta;
  }
}


void NBodySimulation::openParaviewVideoFile () {
  videoFile.open("paraview-output/result.pvd");
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\""
    " version=\"0.1\""
    " byte_order=\"LittleEndian\""
    " compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}

void NBodySimulation::closeParaviewVideoFile () {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
  videoFile.close();
}

void NBodySimulation::printParaviewSnapshot () {
  static int counter = -1;
  counter++;
  std::stringstream filename, filename_nofolder;
  filename << "paraview-output/result-" << counter <<  ".vtp";
  filename_nofolder << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\""
    " NumberOfComponents=\"3\""
    " format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  out.close();

  videoFile << "<DataSet timestep=\"" << counter
            << "\" group=\"\" part=\"0\" file=\"" << filename_nofolder.str()
            << "\"/>" << std::endl;
}

void NBodySimulation::printSnapshotSummary () {
  std::cout << "plot next snapshot"
            << ",\t time step=" << timeStepCounter
            << ",\t t="         << t
            << ",\t dt="        << timeStepSize
            << ",\t v_max="     << maxV
            << ",\t dx_min="    << minDx
            << std::endl;
}

void NBodySimulation::printSummary () {
  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: "
            << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;
}


double* NBodySimulation::remove_element(double arr[], int index, int size) {
    int newSize = size - 1;

    double* newArr = new double[newSize];

    for (int i = 0; i < index; i++) {
        newArr[i] = arr[i];
    }

    for (int i = index + 1; i < size; i++) {
        newArr[i - 1] = arr[i];
    }

    delete[] arr;

    return newArr;
}

double** NBodySimulation::remove_element_2d(double** arr, int index, int size) {

    int newSize = size - 1;

    double** newArr = new double*[newSize];
    for (int i = 0; i < newSize; i++) {
        newArr[i] = new double[3];
    }


    for (int i = 0; i < index; i++) {
        for (int j = 0; j < 3; j++) {
            newArr[i][j] = arr[i][j];
        }
    }


    for (int i = index + 1; i < size; i++) {
        for (int j = 0; j < 3; j++) {
            newArr[i - 1][j] = arr[i][j];
        }
    }


    for (int i = 0; i < size; i++) {
        delete[] arr[i];
    }
    delete[] arr;


    return newArr;
}



double NBodySimulation::calculate_energy() {
    int N = NumberOfBodies;
    double kinetic_energy = 0.0;
    double potential_energy = 0.0;
    for (int i = 0; i < N; i++) {
        kinetic_energy += 0.5 * mass[i] * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
        for (int j = i + 1; j < N; j++) {
            double r = std::sqrt((x[i][0] - x[j][0]) * (x[i][0] - x[j][0]) +
                                 (x[i][1] - x[j][1]) * (x[i][1] - x[j][1]) +
                                 (x[i][2] - x[j][2]) * (x[i][2] - x[j][2]));
            potential_energy -=  mass[i] * mass[j] / r;
        }
    }
    return kinetic_energy + potential_energy;

}

void NBodySimulation::calculate_momentum() {
  double total_momentum[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < NumberOfBodies; i++) {
      total_momentum[0] += mass[i] * v[i][0];
      total_momentum[1] += mass[i] * v[i][1];
      total_momentum[2] += mass[i] * v[i][2];
  }

  std::cout << "Total momentum: (" << total_momentum[0] + total_momentum[1] + total_momentum[2] << ")" << std::endl;
}
