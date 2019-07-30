pipeline {
    agent any

    stages {
        stage('Build') {
            steps {
              scons equations=anisotropic compileMode=release order=6 parallelization=hybrid arch=dsnb compiler=gcc unitTests=fast -j2                
              echo 'Building..'
            }
        }
        stage('Test') {
            steps {
                echo 'Testing..'
            }
        }
        stage('Deploy') {
            steps {
                echo 'Deploying....'
            }
        }
    }
}
