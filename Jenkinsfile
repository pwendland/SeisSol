pipeline {
    agent any

    stages {
        stage('git') {
            steps {
                sh 'git submodule update --init --recursive'
            }
        }
        stage('Build') {
            steps {
              echo 'Building..'
              sh 'which gcc'
              sh 'which mpicc'
              sh 'scons equations=anisotropic compileMode=release order=6 parallelization=hybrid arch=dsnb compiler=gcc -j2'
            }
        }
    }
}
