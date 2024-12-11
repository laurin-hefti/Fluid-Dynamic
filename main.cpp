#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <chrono>

#include "utils.cpp"
#include "grid.cpp"
#include "matrixSolvers.cpp"

class Analyser {
public:
    std::vector<std::vector<float>> saveLog;
    std::chrono::time_point<std::chrono::system_clock> start, end;

    void startTimer(){
        start = std::chrono::system_clock::now();
    }

    void endTimer(){
        end = std::chrono::system_clock::now();
    }

    double getTime(){
        return std::chrono::duration_cast<std::chrono::milliseconds>(this->end - this->start).count();
    }

    void save(std::vector<float> data){
        saveLog.push_back(data);
    }

    void saveToFile(std::string name){
        std::ofstream file(name);

        for (std::vector<float> v : this->saveLog){
            for (float f : v){
                file << f << ",";
            }
            file << ";";
        }

        file.close();
    }
};

Analyser a;

class CFDSolver {
public:
    float tolleranz;
    int maxIter;
    int minIter;

    void init(){
        this->tolleranz = 10e-3;
        this->maxIter = 10e4;
        this->minIter = 10e2;
    }

    std::vector<float> simpleStep(Grid& g, float dt, float fp, float my, int& status){
        //claculate the first velocity
        MatrixSolver solver;

        solver.init();

        g.calc1Ddx2V();
        g.calc1Ddxp();

        std::vector<float> a = g.getMatrixa1D(dt, fp, my);
        std::vector<float> static_a = g.getStaticMatrixa1D();
        std::vector<float> b = g.getMatrixb1D(dt, fp);
        std::vector<float> f_temp(g.numofCells, 0);
        std::vector<float> divergenzV(g.numofCells, 0);
        std::vector<float> p_corr(g.numofCells, 0);
        std::vector<float> p_grad(g.numofCells,0);
        std::vector<float> new_v(g.numofCells, 0);
        std::vector<float> new_p(g.numofCells, 0);

        std::vector<float> old_v(g.numofCells, 0);
        float tolleranz = this->tolleranz;
        float differenz = 1;
        int it = 0;
        
        while(differenz > tolleranz || it < this->minIter){
            
            g.calc1Ddx2V();
            g.calc1Ddxp();
            
            a = g.getMatrixa1D(dt, fp, my);
            //a = g.getStaticMatrixa1D();
            
            b = g.getMatrixb1D(dt, fp);

            f_temp = solver.GausSeidel(a,b);
            
            //calculate the preasure correcture
            
            scaleVec(f_temp,fp/dt);
            
            //g.set1Dvel(f_temp); //?
            
            g.calc1DVeldiv();
            
            divergenzV = g.getVelDivergenz();
            
            p_corr = solver.Jakobi(static_a,divergenzV);
            
            //add the preasure corr to the vell
            scaleVec(p_corr, 0.6);
            
            g.addPreasure(p_corr);
            
            //calc the p gradient
            g.calc1Ddxp();
            
            p_grad = g.getPresGrad();
            
            scaleVec(p_grad, dt/fp);
            
            //sub the p gradient to the vel
            
            
            //scaleVec(p_grad,0.8);
            
            g.subVelocity(p_grad);
            
            //calc div
            
            new_v = g.getVel();
            new_p = g.getPres();
            
            differenz = 0;
            
            for (int i = 0; i < new_v.size(); i++){
                if (differenz < std::abs(new_v[i]-old_v[i])){
                    differenz = std::abs(new_v[i]-old_v[i]);
                }
                
                old_v[i] = new_v[i];
            }
            
            //g.set1Dvel(old_v); //not needed because vel iun the grid willb e corrected with the correcture
            
            it += 1;
            if (it > this->maxIter){
                std::cout << "not conv simple cfd solver dif : " << differenz << std::endl;
                status = false;
                break;
            }
        }
        
        return old_v;
    }
    
    void simplerun(Grid& g, float startT, float endT, float dt, float fp, float my){
        
        int status = true;

        std::vector<float> adv_pres(g.numofCells, 0);
        std::vector<float> adv_vel(g.numofCells, 0);
        std::vector<float> new_v(g.numofCells, 0);
        
        while (status && startT < endT){
            
            new_v = this->simpleStep(g, dt, fp, my, status);

            g.set1Dvel(new_v); // should be deminsionsles
            
            //std::vector<float> new_p = g.getPres();

            g.calc1Ddxp();
            g.calc1Ddxv();

            adv_pres = g.calc1DpresAdv(dt);
            adv_vel = g.calc1DvelAdv(dt);

            g.set1Dpres(adv_pres);
            g.set1Dvel(adv_vel);

            //std::vector<float> veldiv = g.getVelDivergenz();

            a.save(adv_pres);

            g.zeroAllGrad1D();
            
            startT += dt;
        }
        
        if (!status){
            std::cout << "error at time : " << startT << std::endl;
        } else {
            std::cout << "finish simulation" << std::endl;
        }
        
        printMatrix(g.getPres(),20);
    }
};

int main() {

    a.startTimer();
    Grid g;
    g.init1D(100, 1,1);

    //std::vector<float> initPfield = {1,1,1,0.75,0.5,0.25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    std::vector<float> initPfield = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.3,0.5,0.8,1,0.8,0.5,0.3,0.1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    std::vector<float> initFfield = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    g.set1Dpres(initPfield);
    g.set1Dvel(initFfield);
    
    CFDSolver sys;
    sys.init();

    sys.simplerun(g, 0,4,0.01, 1000, 10e-4);

    a.endTimer();

    std::cout << "time : " << a.getTime();
    
    a.saveToFile("log.txt");
    return 0;
}
