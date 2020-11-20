#ifndef POTENTIAL_H
#define POTENTIAL_h
class Sim;

class Potential{
    private:
        double omega;
        Sim *sim;
        double envelope(double x, double z, double t);
        
    public:
        Potential(Sim *_sim);
        double potential(double x, double z, double t);
};

#endif