class MatrixSolver {
public:
    std::vector<float> Jakobi(std::vector<float>& a, std::vector<float>& b){
        int xsize = b.size();
        int ysize = a.size() / xsize;
        
        float tolleranz = 10e-4;
        float diverenz = 1;
        int maxiter = 10e5;
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
                
                if (diverenz > 10e7 || it < 1000 || dump){
                    x[i] = x[i]/2 + new_val/2;
                    
                    if (diverenz < 10e4){dump = false;}
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
            
            if (diverenz > 10e9){
                std::cout << "jakobi solver error to big" << std::endl;
                break;
            }
            
            it += 1;

            if (it > maxiter){
                std::cout << "jakobi solver not converged error : " << diverenz << std::endl;
                break;
            }
        }
        
        return x;
    }
    
};