#define TEST false

#if TEST == true
#include <iostream>
#include <vector>
#include "utils.cpp"
#endif

class MatrixSolver {
public:
float tolleranz;
int maxiter;
int startDump;
float dumpinit;
float dumpreleas;
float maxError;

    void init(){
        this->tolleranz = 10e-4;
        this->maxiter = 10e5;
        this->startDump = 1000;
        this->dumpinit = 10e7;
        this->dumpreleas = 10e4;
        this->maxError = 10e9;
    }

    std::vector<float> Jakobi(const std::vector<float>& a, const std::vector<float>& b){
        int xsize = b.size();
        int ysize = a.size() / xsize;
        
        float tolleranz = this->tolleranz;
        float diverenz = 1;
        int maxiter = this->maxiter;
        std::vector<float> oldx(xsize, 0);
        std::vector<float> x(xsize,0);
        
        bool dump = false;
        
        int it = 0;
        
        while (diverenz > tolleranz || it < 1) {
            for (int i = 0; i < xsize; i++){
                float new_val = 1 / a[i*ysize + i];
                float temp_val = b[i];
                
                for (int j = 0; j < ysize; j++){
                    if (i != j) {
                        temp_val -= a[i*ysize+j] * oldx[j]; //mustz be oldx[j]
                    }
                }
                
                new_val *= temp_val;
                
                if (diverenz > this->dumpinit || it < this->startDump || dump){
                    x[i] = x[i]/2 + new_val/2;
                    
                    if (diverenz < this->dumpreleas){dump = false;}
                } else {
                    x[i] = new_val;
                }
            }
            
            diverenz = 0;
            for(int z = 0; z < xsize; z++){
                float t = std::abs(oldx[z]-x[z]);
                if (t > diverenz){
                    diverenz = t;
                }
                
                oldx[z] = x[z];
            }
            
            if (diverenz > this->maxError){
                std::cout << "jakobi solver error to big" << std::endl;
                break;
            }
            
            it += 1;

            if (it > maxiter){
                std::cout << "jakobi solver not converged error : " << diverenz << std::endl;
                break;
            }
        }

        #if TEST == true
        std::cout << "jakob iteration : " << it << std::endl;
        #endif

        return x;
    }

    std::vector<float> GausSeidel(std::vector<float>& a, std::vector<float>& b){
        int xsize = b.size();
        int ysize = a.size() / xsize;
        
        float tolleranz = this->tolleranz;
        float diverenz = 1;
        int maxiter = this->maxiter;
        std::vector<float> oldx(xsize, 0);
        std::vector<float> x(xsize,0);
        
        bool dump = false;
        
        int it = 0;
        
        while (diverenz > tolleranz || it < 1) {
            for (int i = 0; i < xsize; i++){
                float new_val = 1 / a[i*ysize + i];
                float temp_val = b[i];
                
                for (int j = 0; j < i; j++){
                    if (i != j) {
                        temp_val -= a[i*ysize+j] * x[j]; //mustz be oldx[j]
                    }
                }

                for (int j = i; j < ysize; j++){
                    if (i != j) {
                        temp_val -= a[i*ysize+j] * oldx[j]; //mustz be oldx[j]
                    }
                }
                
                new_val *= temp_val;
                
                if (diverenz > this->dumpinit || it < this->startDump || dump){
                    x[i] = x[i]/2 + new_val/2;
                    
                    if (diverenz < this->dumpreleas){dump = false;}
                } else {
                    x[i] = new_val;
                }
            }

            diverenz = 0;
            for(int z = 0; z < xsize; z++){
                float t = std::abs(oldx[z]-x[z]);
                if (t > diverenz){
                    diverenz = t;
                }
                
                oldx[z] = x[z];
            }
            
            if (diverenz > this->maxError){
                std::cout << "jakobi solver error to big" << std::endl;
                break;
            }
            
            it += 1;

            if (it > maxiter){
                std::cout << "jakobi solver not converged error : " << diverenz << std::endl;
                break;
            }
        }
        
        #if TEST == true
        std::cout << "Gaus Seidel iteration : " << it << std::endl;
        #endif
        
        return x;
    }
    
};

#if TEST == true
int main(){
    MatrixSolver solver;

    solver.init();

    std::vector<float> a = {2,-1,0,-1,2,-1,0,-1,2};
    std::vector<float> b = {0.1,0.005,0.1};

    std::vector<float> x_jakobi = solver.Jakobi(a,b);

    std::vector<float> x_gausSeidel = solver.GausSeidel(a,b);
    printMatrix(x_jakobi,3);
    printMatrix(x_gausSeidel,3);
    return 0;
}
#endif