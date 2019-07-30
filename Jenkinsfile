pipeline {
    agent any

    stages {
        stage('Build') {
            steps {
              echo 'Building..'
              scons equations=elastic compileMode=release order=6 parallelization=hybrid arch=dsnb compiler=gcc unitTests=fast -j2                
            }
        }
    }
}
