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

double NBodySimulation::force_calculation (int i, int j, int direction){
  // Euclidean distance
  const double distance = sqrt(
                               (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
                               (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
                               (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
                               );;
  const double distance3 = distance * distance * distance;
  minDx = std::min( minDx,distance );

  return (x[i][direction]-x[j][direction]) * mass[i]*mass[j] / distance3;
}

double NBodySimulation::NBodySimulationMolecularForces(int i, int j, int direction) {
    double rij[3]; 
    double rij_norm; 
    double force_ij; 
    double direction_vector; 

 
    for (int d = 0; d < 3; d++) {
        rij[d] = x[i][d] - x[j][d];
    }

    rij_norm = sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
    

    force_ij = 10 * ((0.1 / pow(rij_norm, 13)) - (0.1 / pow(rij_norm, 9))) * rij_norm;

    direction_vector = force_ij * rij[direction] / rij_norm;

    return direction_vector;
}

void NBodySimulation::updateBody () {
  // std::cout << "energy: "<< calculate_energy() <<std::endl;
  timeStepCounter++;
  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();

  double half_dt = timeStepSize/2;

// compute half an Euler time step for dV --------------------------------------
  double* force0 = new double[NumberOfBodies]{0,0};
  double* force1 = new double[NumberOfBodies]{0,0};
  double* force2 = new double[NumberOfBodies]{0,0};

    for (int i = 0; i < NumberOfBodies; i++) {
      for (int j = i+1; j < NumberOfBodies; j++) {

          double Fx = force_calculation(j, i, 0);
          double Fy = force_calculation(j, i, 1);
          double Fz = force_calculation(j, i, 2);


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

          double Fx = force_calculation(j, i, 0);
          double Fy = force_calculation(j, i, 1);
          double Fz = force_calculation(j, i, 2);


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
   double C = 0.01/NumberOfBodies;
    for (int i = 0; i < NumberOfBodies; i++) {
      for (int j = i+1; j < NumberOfBodies; j++) {
          // std::cout <<euclidean_distance(i,j) <<"   " <<C*(mass[i]+mass[j]) <<std::endl;
          // std::cout <<i<<"   " <<j <<std::endl;
          if (euclidean_distance(i,j) <= C*(mass[i]+mass[j])){
            std::cout << "colled " <<std::endl;
            //do coll
            for (int direction =0;direction <= 2 ;direction++ ){
              x[i][direction]= (mass[i]*x[i][direction]+mass[j]*x[j][direction])/(mass[i]+mass[j]);
              v[i][direction]= (mass[i]*v[i][direction]+mass[j]*v[j][direction])/(mass[i]+mass[j]);
              
            }
            mass[i] = mass[i] + mass[j];
            x[j] = x[NumberOfBodies-1];
            v[j] = v[NumberOfBodies-1];
            mass[j] = mass[NumberOfBodies-1];
            NumberOfBodies--;
            
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
