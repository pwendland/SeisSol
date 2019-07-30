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
              sh 'scons equations=elastic compileMode=release order=6 parallelization=hybrid netcdf=yes arch=dsnb compiler=gcc -j2'
              sh 'cp build/SeisSol_release_generatedKernels_dsnb_hybrid_none_9_6 seissol-benchmarks/SeisSol'
              sh 'echo $PWD/Maple/ > seissol-benchmarks/DGPATH'
            }
        }
        stage('Run') {
            steps {
              echo 'Run..'
              sh 'cd seissol-benchmarks && ./SeisSol loh1_parameters.par'
              sh 'diff seissol-benchmarks/output/LOH1-receiver-00001-00000.dat seissol-benchmarks/output/reference/coarse/LOH1-receiver-00001-00000.dat'
            }
        }
    }
}
