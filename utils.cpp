void scaleVec(std::vector<float>& v, float s){
    for (int i = 0; i < v.size(); i++){
        v[i] *= s;
    }
}

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