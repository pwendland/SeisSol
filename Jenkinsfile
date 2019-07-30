pipeline {
    agent any

    stages {
        stage('git') {
            steps {
              sh 'git submodule update --init --recursive'
              sh 'git clone git@gitlab.lrz.de:ga35kum/seissol-benchmarks.git'
            }
        }
        stage('Build') {
            steps {
              echo 'Building..'
              sh 'which gcc'
              sh 'which mpicc'
              sh 'scons equations=elastic compileMode=release order=6 parallelization=hybrid arch=dsnb compiler=gcc -j2'
              sh 'cp build/build_release_generatedKernels_dsnb_hybrid_none_9_6 seissol-benchmarks/SeisSol'
              sh 'echo $PWD/Maple/ > seissol-benchmarks/DGPATH'
            }
        }
        stage('Run') {
            echo 'Run..'
            sh 'cd seissol-benchmarks'
            sh './SeisSol loh1_parameters.par'
        }
    }
}
