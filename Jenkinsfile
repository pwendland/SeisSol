pipeline {
    agent any

    stages {
        stage('Build') {
            steps {
              echo 'Building..'
              scons equations=anisotropic compileMode=release order=6 parallelization=hybrid arch=dsnb compiler=gcc -j2
            }
        }
    }
}
