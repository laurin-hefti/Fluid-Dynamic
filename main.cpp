#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

void scaleVec(std::vector<float>& v, float s){
    for (int i = 0; i < v.size(); i++){
        v[i] *= s;
    }
}

void printMatrix(std::vector<float> m, int x);

typedef struct Cell {
    float vel;
    int vel_dir;
    float pres;
    
    float veldiv; //velocity divergenz
    
    float dxp;
    float dxv;
    float dx2v;
} Cell;

class Grid {
public:
    Cell* cells;
    int numofCells;
    int d;
    
    float dx;
    
    int border_correcture; //0 = none, 1 = neuman
                            //2 = dietrich // periodic
    
    void init1D(int xs, float dx){
        this->d = 1;
        this->numofCells = xs;
        this->dx = dx;
        
        this->border_correcture = 1;
        
        this->cells = (Cell*) malloc(sizeof(Cell) * (xs+1));
        
        for (int i = 0; i < xs; i++){
            this->cells[i] = Cell{0,0,0,0,0,0};
        }
        
        this->cells[xs] = Cell{0,0,0,0,0,0};
        
    }
    
    int corr(int i){
        if (this->border_correcture == 0) {return i;}
        else if (this->border_correcture == 1){ //neuman corr
            if (i < 0) {return i+1;}
            else if (i >= this->numofCells) {return this->numofCells-1;}
            else {return i;}
        }
        else if (this->border_correcture == 2){
            if (i < 0 || i >= this->numofCells){
                return this->numofCells;
            } else {
                return i;
            }
        }
        else if (this->border_correcture == 3){
            if (i < 0) {return this->numofCells-1;}
            else if (i >= this->numofCells) {return 0;}
            else {return i;}
        }
        else {return -1;}
    }
    
    void set1Dpres(std::vector<float> p){
        for (int i = 0; i < this->numofCells; i++){
            this->cells[i].pres = p[i];
        }
    }
    
    
    void set1Dvel(std::vector<float> v){
        for (int i = 0; i < this->numofCells; i++){
            this->cells[i].vel = v[i];
        }
    }
    
     void calc1Ddx2V(){
         for (int i = 0; i < numofCells; i++){
             float xi_1 = this->cells[this->corr(i-1)].vel;
             float xi = this->cells[this->corr(i)].vel;
             float xi1 = this->cells[this->corr(i+1)].vel;
             
             this->cells[i].dx2v = (xi_1 + xi1 - 2* xi) / (this->dx * this->dx);
         }
     }
     
     void calc1Ddxp(){
         for (int i = 0; i < numofCells; i++){
             float x1 = this->cells[this->corr(i+1)].pres;
             float x_1 = this->cells[this->corr(i-1)].pres;
             
             this->cells[i].dxp = (x1 - x_1) / (2 * this->dx);// * this->dx);
         }
     }
     
     void calc1Ddxv(){
         for (int i = 0; i < numofCells; i++){
             float x1 = this->cells[this->corr(i+1)].vel;
             float x_1 = this->cells[this->corr(i-1)].vel;
             
             this->cells[i].dxv = (x1 - x_1) / (2 * this->dx);// * this->dx);
         }
     }
     
     void calc1DVeldiv(){
         for (int i = 0; i < numofCells; i++){
             float x1 = this->cells[this->corr(i+1)].vel;
             float x_1 = this->cells[this->corr(i-1)].vel;
             
             this->cells[i].veldiv = (x1 - x_1) / (2 * this->dx);
         }
     }
     
     void addPreasure(std::vector<float> p){
         for (int i = 0; i < this->numofCells; i++){
             this->cells[i].pres += p[i];
         }
     }
     
     void subVelocity(std::vector<float> v){
         for (int i = 0; i < this->numofCells; i++){
             this->cells[i].vel += v[i];                //eigentlich -= v[i];
         }
     }
     
     std::vector<float> getPresGrad(){
         std::vector<float> pg;
         for (int i = 0; i < this->numofCells; i++){
             pg.push_back(this->cells[i].dxp);
         }
         
         return pg;
    }
    
    std::vector<float> getPres(){
         std::vector<float> p;
         for (int i = 0; i < this->numofCells; i++){
             p.push_back(this->cells[i].pres);
         }
         
         return p;
    }
    
     std::vector<float> getVel(){
         std::vector<float> v;
         for (int i = 0; i < this->numofCells; i++){
             v.push_back(this->cells[i].vel);
         }
         
         return v;
    }
    
     std::vector<float> getVelDivergenz(){
         std::vector<float> v;
         for (int i = 0; i < this->numofCells; i++){
             v.push_back(this->cells[i].veldiv);
         }
         
         return v;
    }
     
     std::vector<float> getMatrixb1D(float dt, float fp){
         std::vector<float> b;
         
         for (int i = 0; i < this->numofCells; i++){
             b.push_back(
                 (this->cells[i].vel / dt) - 
                 (1/fp * ((this->cells[this->corr(i+1)].pres - this->cells[this->corr(i-1)].pres) / this->dx))// + this->cells[i].vel
             ); //external forces
         }
         
         return b;
    }
     
     std::vector<float> getMatrixa1D(float dt, float fp, float my){ //very standart
        std::vector<float> a;
        for (int i = 0; i < this->numofCells; i++){
            for (int j = 0; j < this->numofCells; j++){
                if (i == j){
                    a.push_back((fp/dt) + (2*this->cells[i].vel)/(this->dx*this->dx));
                } else if (j+1 == i || j-1 == i){
                    a.push_back(-(my/fp) / (this->dx*this->dx));
                }
                else {
                    a.push_back(0);
                }
            }
        }
     return a;
     }
};

class MatrixSolver {
public:
    std::vector<float> Jakobi(std::vector<float> a, std::vector<float> b){
        int xsize = b.size();
        int ysize = a.size() / xsize;
        
        float tolleranz = 10e-5;
        float diverenz = 1;
        int maxiter = 10e6;
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
                    
                    if (diverenz < 10e2){dump = false;}
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

class CFDSolver {
public:
    std::vector<float> simpleStep(Grid& g, float dt, float fp, float my, int& status){
        //claculate the first velocity
        MatrixSolver solver;
        
        std::vector<float> a = g.getMatrixa1D(dt, fp, my);
        
        std::vector<float> old_v(g.numofCells, 0);
        float tolleranz = 0.01;
        float differenz = 1;
        int it = 0;
        
        while(differenz > tolleranz || it < 100){
            
            g.calc1Ddx2V();
            g.calc1Ddxp();
            
            a = g.getMatrixa1D(dt, fp, my);
            
            std::vector<float> b = g.getMatrixb1D(dt, fp);
            //printMatrix(a,6);
            std::vector<float> f_temp = solver.Jakobi(a,b);
            
            //calculate the preasure correcture
            
            scaleVec(f_temp,fp/dt);
            
            //g.set1Dvel(f_temp); //?
            
            g.calc1DVeldiv();
            
            std::vector<float> divergenzV = g.getVelDivergenz();
            
            std::vector<float> p_corr = solver.Jakobi(a,divergenzV);
            
            //add the preasure corr to the vell
            scaleVec(p_corr, 0.8);
            
            g.addPreasure(p_corr);
            
            //calc the p gradient
            g.calc1Ddxp();
            
            std::vector<float> p_grad = g.getPresGrad();
            
            scaleVec(p_grad, dt/fp);
            
            //printMatrix(p_grad,4);
            
            //sub the p gradient to the vel
            
            
            //scaleVec(p_grad,0.8);
            
            g.subVelocity(p_grad);
            
            //calc div
            
            std::vector<float> new_v = g.getVel();
            std::vector<float> new_p = g.getPres();
            
            if (it > 1000){
            //printMatrix(new_p,6);
            }
            
            differenz = 0;
            
            for (int i = 0; i < new_v.size(); i++){
                if (differenz < std::abs(new_v[i]-old_v[i])){
                    differenz = std::abs(new_v[i]-old_v[i]);
                }
                
                old_v[i] = new_v[i];
            }
            
            //g.set1Dvel(old_v); not needed because vel iun the grid willb e corrected with the correcture
            
            it += 1;
            if (it > 6000){
                std::cout << "not conv simple cfd solver dif : " << differenz << std::endl;
                status = false;
                break;
            }
            //end loop
        }
        
        return old_v;
    }
    
    void simplerun(Grid& g, float startT, float endT, float dt){
        float fp = 1000; //fluid preasure
        float my = 10e-3;
        
        int status = true;
        
        while (status && startT < endT){
            
            std::vector<float> new_v = this->simpleStep(g, dt, fp, my, status);
            
            g.set1Dvel(new_v); // should be deminsionsles
            
            std::vector<float> new_p = g.getPres();
            
            //printMatrix(new_p,6);
            startT += dt;
        }
        
        if (!status){
            std::cout << "error at time : " << startT << std::endl;
        } else {
            std::cout << "finish simulation" << std::endl;
        }
        
        printMatrix(g.getPres(),6);
    }
};

void printMatrix(std::vector<float> m, int x){
    int len = m.size();
    int y = len/x;
    for (int i = 0; i < y; i++){
        for (int j = 0; j < x; j++){
            std::cout << m[i*x+j] << " ";
        }
        std::cout << std::endl;
    }
}


int main() {
    Grid g;
    g.init1D(6, 0.5);
    
    std::vector<float> initPfield = {1,0.5,0,0,0,0};
    std::vector<float> initFfield = {0,0,0,0,0,0};
    
    g.set1Dpres(initPfield);
    g.set1Dvel(initFfield);
    
    CFDSolver sys;
    
    sys.simplerun(g, 0,40,0.1);
    //int status = true;
    // grid, dt, fp, my
    //std::vector<float> new_v = sys.simpleStep(g, 0.1, 1000, 10e-3, status);
    //std::cout << "solved : " << std::endl;
    //printMatrix(new_v, 6);
    /*
    MatrixSolver solver;
    
    std::vector<float> a = {4,1,2,1,3,2,1,1,2};
    std::vector<float> b = {12,13,9};
    
    std::vector<float> x = solver.Jakobi(a,b);
    
    printMatrix(x, 3);
    */
    
    return 0;
}
