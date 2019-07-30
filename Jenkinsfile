pipeline {
    agent any

    stages {
        stage('Build') {
            steps {
              echo 'Building..'
              which gcc
              which mpicc
              scons equations=anisotropic compileMode=release order=6 parallelization=hybrid arch=dsnb compiler=gcc -j2
            }
        }
    }
}
