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
    
    int border_correcture;

    /* 
    0 = none, 1 = dv/dx = 0, neumann, 2 = perodic 
    */
    
    void init1D(int xs, float dx, int border_correcture){
        this->d = 1;
        this->numofCells = xs;
        this->dx = dx;
        
        this->border_correcture = border_correcture;
        
        this->cells = (Cell*) malloc(sizeof(Cell) * (xs+2));
        
        for (int i = 0; i < xs; i++){
            this->cells[i] = Cell{0,0,0,0,0,0,0};
        }
        
        this->cells[xs] = Cell{0,0,0,0,0,0,0};
        this->cells[xs+1] = Cell{0,0,0,0,0,0,0};
        
    }

    /*Cell& requestCell(int i){

    }*/
    
    int corr(int i){
        if (this->border_correcture == 0) {return i;}

        else if (this->border_correcture == 1){
            if (i >= this->numofCells) {return this->numofCells+1;}
            else if (i < 0) {return this->numofCells;}
            else {return i;}
        }

        else if (this->border_correcture == 2){
            if (i < 0) {return this->numofCells-1;}
            else if (i >= this->numofCells) {return 0;}
            else {return i;}
        }

        std::cout << "error not border correction" << std::endl;
        return -1;
    }

    void copyCell(int i, int j){
        this->cells[i].vel = this->cells[j].vel;
        this->cells[i].pres = this->cells[j].pres;
        this->cells[i].dxp = this->cells[j].dxp;
        this->cells[i].dxv = this->cells[j].dxv;
        this->cells[i].dx2v = this->cells[j].dx2v;
    }

    void updateBorder(){
        if (this->border_correcture == 1){
            this->copyCell(this->numofCells, 0);
            this->copyCell(this->numofCells+1, this->numofCells-1);
            this->cells[this->numofCells].vel = 0;
            this->cells[this->numofCells+1].vel = 0;    
        }
    }
    
    void set1Dpres(const std::vector<float>& p){
        for (int i = 0; i < this->numofCells; i++){
            this->cells[i].pres = p[i];
        }

        this->updateBorder();
    }
    
    
    void set1Dvel(const std::vector<float>& v){
        for (int i = 0; i < this->numofCells; i++){
            this->cells[i].vel = v[i];
        }

        this->updateBorder();
    }
    
     void calc1Ddx2V(){
         for (int i = 0; i < numofCells; i++){
             float xi_1 = this->cells[this->corr(i-1)].vel;
             float xi = this->cells[this->corr(i)].vel;
             float xi1 = this->cells[this->corr(i+1)].vel;
             
             this->cells[i].dx2v = (xi_1 + xi1 - 2* xi) / (this->dx * this->dx);
         }

         this->updateBorder();
     }
     
     void calc1Ddxp(){
         for (int i = 0; i < numofCells; i++){
             float x1 = this->cells[this->corr(i+1)].pres;
             float x_1 = this->cells[this->corr(i-1)].pres;
             
             this->cells[i].dxp = (x1 - x_1) / (2 * this->dx);// * this->dx);
         }

         this->updateBorder();
     }
     
     void calc1Ddxv(){
         for (int i = 0; i < numofCells; i++){
             float x1 = this->cells[this->corr(i+1)].vel;
             float x_1 = this->cells[this->corr(i-1)].vel;
             
             this->cells[i].dxv = (x1 - x_1) / (2 * this->dx);// * this->dx);
         }

         this->updateBorder();
     }
     
     void calc1DVeldiv(){
         for (int i = 0; i < numofCells; i++){
             float x1 = this->cells[this->corr(i+1)].vel;
             float x_1 = this->cells[this->corr(i-1)].vel;
             
             this->cells[i].veldiv = (x1 - x_1) / (2 * this->dx);
         }

         this->updateBorder();
     }
     
     void addPreasure(const std::vector<float>& p){
         for (int i = 0; i < this->numofCells; i++){
             this->cells[i].pres += p[i];
         }

         this->updateBorder();
     }
     
     void subVelocity(const std::vector<float>& v){
         for (int i = 0; i < this->numofCells; i++){
             this->cells[i].vel -= v[i];                //eigentlich -= v[i];
         }

         this->updateBorder();
     }
     
     std::vector<float> getPresGrad(){
         std::vector<float> pg(this->numofCells, 0);
         for (int i = 0; i < this->numofCells; i++){
             pg[i] = (this->cells[i].dxp);
         }
         
         return pg;
    }
    
    std::vector<float> getPres(){
         std::vector<float> p(this->numofCells, 0);
         for (int i = 0; i < this->numofCells; i++){
             p[i] = (this->cells[i].pres);
         }
         
         return p;
    }
    
     std::vector<float> getVel(){
         std::vector<float> v(this->numofCells, 0);
         for (int i = 0; i < this->numofCells; i++){
             v[i] = (this->cells[i].vel);
         }
         
         return v;
    }
    
     std::vector<float> getVelDivergenz(){
         std::vector<float> v(this->numofCells, 0);
         for (int i = 0; i < this->numofCells; i++){
             v[i] = (this->cells[i].veldiv);
         }
         
         return v;
    }
     
     std::vector<float> getMatrixb1D(float dt, float fp){
         std::vector<float> b(this->numofCells, 0);
         
         for (int i = 0; i < this->numofCells; i++){
             b[i] = (
                 (this->cells[i].vel / dt) - 
                 (1/fp * 
                 ((this->cells[this->corr(i+1)].pres - this->cells[this->corr(i-1)].pres) 
                 / (this->dx*2)))// + this->cells[i].vel //noral nicht mal zwei
             ); //external forces
         }
         
         return b;
    }
     
     std::vector<float> getMatrixa1D(float dt, float fp, float my){ //very standart
        std::vector<float> a(this->numofCells * this->numofCells, 0);
        for (int i = 0; i < this->numofCells; i++){
            for (int j = 0; j < this->numofCells; j++){
                if (i == j){
                    a[i*this->numofCells + j] = ((fp/dt) + 
                    (2*this->cells[i].vel)
                    /(this->dx*this->dx));
                } else if (j+1 == i || j-1 == i){
                    a[i*this->numofCells + j] = (-(my/fp) / (this->dx*this->dx));
                }
                else {
                    a[i*this->numofCells + j] = 0;
                }
            }
        }
     return a;
     }

     std::vector<float> getStaticMatrixa1D(){
        std::vector<float> a(this->numofCells * this->numofCells, 0);
        for (int i = 0; i < this->numofCells; i++){
            for (int j = 0; j < this->numofCells; j++){
                if (i == j){
                    a[i*this->numofCells + j] = -2;
                } else if (j+1 == i || j-1 == i){
                    a[i*this->numofCells + j] = 1;
                }
                else {
                    a[i*this->numofCells + j] = 0;
                }
            }
        }
     return a;
     }

     std::vector<float> calc1DpresAdv(float dt){
        std::vector<float> adv(this->numofCells, 0);

        this->calc1Ddxp();

        for(int i = 0; i < this->numofCells; i++){
            Cell c = this->cells[i];
            adv[i] = (
                c.pres - (dt * c.vel * c.dxp)
            );
        }

        return adv;
     }

     std::vector<float> calc1DvelAdv(float dt){
        std::vector<float> adv(this->numofCells, 0);

        this->calc1Ddxv();

        for (int i = 0; i < this->numofCells; i++){
            Cell c = this->cells[i];
            adv[i] = (
                c.vel - (dt * c.vel * c.dxv) // + external forces
            );
        }

        return adv;
     }

     void zeroAllGrad1D(){
        for (int i = 0; i < this->numofCells; i++){
            this->cells[i].dxp = 0;
            this->cells[i].dxv = 0;
            this->cells[i].dx2v = 0;
        }

        this->updateBorder();
     }
};