#ifndef POTENTIAL_H
#define POTENTIAL_h
class Sim;

class Potential{
    private:
        double omega;
	double z0;
        Sim *sim;
        double envelope(double x, double z);
        
    public:
        Potential(Sim *_sim);
        double potential(double x, double z, double t);
};

#endif
