#include <memory>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <cmath>
#include <map>
#include <numeric>

using namespace std;

class Measurer;

class IsingMCWorker {
public:
    IsingMCWorker(int _Nx, int _Ny, double _J);

    void init_lattice_random();
    void print_system();
    
    inline short& lat(int x, int y) {
        // Periodic boundary conditions
        return lattice[((x + Nx)%Nx)*Ny + (y+Ny)%Ny];
    }

    void metropolis(int Nsteps, double T, Measurer *m);
    // Lattice parameters
    int Nx, Ny, Nsites;
    vector<short> lattice;
    // System parameters
    double J;
    // Observables
    double energy;
    int magnetization;
};

class Measurer {
public:
    Measurer(IsingMCWorker *_system);
    void start();
    void measure();

    void calc_results();
    void print_results();

    void spin_flip(int x, int y, int time);
    double calc_autocorr(int maxtime, int diff);

    static int dict_sum_d(vector<int> times, int maxtime, int diff);
private:
    IsingMCWorker* system;
    bool active = false;
    vector<short> init_lattice;

    vector<int> magnetizations;
    vector<double> energies;

    double energy = 0.0;
    double magnetization = 0.0;

    vector<vector<int>> flip_times;
};

IsingMCWorker::IsingMCWorker(int _Nx, int _Ny, double _J)
    : Nx(_Nx), Ny(_Ny), Nsites(Nx*Ny), J(_J)
{
    energy = 0.0;
    magnetization = 0.0;
    lattice.resize(Nsites);
}

void IsingMCWorker::init_lattice_random() 
{
    magnetization = 0;
    generate(lattice.begin(), lattice.end(), [&]() {
                short s = 2 * (rand() % 2) - 1;
                magnetization += s;
                return s;
            });

    energy = 0.0;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            energy += -J*lat(i, j)*(lat(i, j+1) + lat(i+1,j));
        }
    }
}

void IsingMCWorker::print_system() {
    cout << " === System state === " << endl;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            cout << (lat(i,j) == +1 ? '+' : '-');
        }
        cout << endl;
    }
    cout << "Magnetization: " << magnetization << endl;
    cout << "Energy: " << energy << endl;
}

void IsingMCWorker::metropolis(int time, double T, Measurer *m = nullptr) 
{
	for (long int i = 0; i < time; i++)
	{
        for (int j = 0; j < Nsites; j++) {
            int x = rand() % Nx;
            int y = rand() % Ny;
            int sum = lat(x-1,y) + lat(x+1, y) + lat(x,y-1) + lat(x, y+1);
            
            double dE = 2*J*lat(x,y)*sum;

            if (dE <= 0 || (rand() / (double)RAND_MAX) < exp(-dE/T))
            {
                lat(x, y) = -lat(x, y);
                magnetization += 2*lat(x,y);
                energy += dE;
                if (m) m->spin_flip(x, y, i*Nsites + j);
            }
        }
        if (m) m->measure();
	}
}

Measurer::Measurer(IsingMCWorker *_system) 
    : system(_system)
{
}

void Measurer::start()
{
    active = true;
    init_lattice = system->lattice;
    flip_times.resize(system->Nsites);
}

void Measurer::measure()
{
    if (!active) return;
    magnetizations.push_back(system->magnetization);
    energies.push_back(system->energy);
}

void Measurer::spin_flip(int x, int y, int time) {
    if (!active) return;
    flip_times[x*system->Ny + y].push_back(time);
}

int Measurer::dict_sum_d(vector<int> times, int maxtime, int diff)
{
    int sum = 0;
    vector<int> timesSort;

    for (int i = 0; i < times.size(); i++) {
        if (times[i] >= diff) {
            timesSort.push_back(times[i]);
        }
        if (times[i] + diff < maxtime) {
            timesSort.push_back(times[i] + diff);
        }
    }

    timesSort.push_back(maxtime);
    sort(timesSort.begin(), timesSort.end());
    
    int curt = diff;
    short flag = +1;
    for (int i = 0; (i < timesSort.size()) && (timesSort[i] < diff); i++)
        flag = -flag;
    for (auto p : timesSort) {
        sum += (p - curt)*flag;
        flag = -flag;
        curt = p;
    }

    return sum;
}

void Measurer::calc_results()
{
    energy = accumulate(energies.begin(), energies.end(), 0.0);
    energy /= energies.size();
    energy /= system->Nsites;

    magnetization = accumulate(magnetizations.begin(), magnetizations.end(), 0.0);
    magnetization /= magnetizations.size();
    magnetization /= system->Nsites;
}

double Measurer::calc_autocorr(int maxtime, int diff)
{
    double sum = 0;
    for (auto v : flip_times) 
    {
        sum += Measurer::dict_sum_d(v, maxtime, diff);
    }
    sum /= (maxtime - diff);
    sum /= system->Nsites;
    return sum - magnetization*magnetization;
}


void Measurer::print_results() 
{
    cout << "Energies: ";
    for (auto E : energies) {
        cout << (double)E / system->Nsites << " ";
    }
    cout << endl;

    cout << "Magnetizations: ";
    for (auto M : magnetizations) {
        cout << (double)M / system->Nsites << " ";
    }
    cout << endl;
}

int main()
{
    double J = 1.0;
    double temp = 1.0;
    int Ttherm = 10000;
    int Trun = 10000;


    srand(time(nullptr));

    IsingMCWorker mc(100, 100, J);
    mc.init_lattice_random();

    cout << "Thermalizing...";
    Measurer m(&mc);
    mc.metropolis(Ttherm, temp, &m);
    cout << "done" << endl;

    mc.print_system();
   
    cout << "Measuring...";
    m.start();
    mc.metropolis(Trun, temp, &m);
    cout << "done" << endl;

    m.calc_results();
    for (int t = 0; t < Trun / 2; t += 10) {
        cout << t << " " << m.calc_autocorr(Trun*mc.Nsites, t*mc.Nsites) << endl;
    }

    return 0;
}
