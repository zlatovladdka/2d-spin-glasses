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
    void measure();

    void calc_results();
    void print_results();

    void spin_flip(int x, int y, int time);
    void calc_avg(int maxtime);
    double calc_autocorr(int maxtime, int diff);

    double energy = 0.0;
    double magnetization = 0.0;
private:
    IsingMCWorker* system;

    vector<int> magnetizations;
    vector<double> energies;

    vector<vector<int>> flip_times;
    vector<double> spin_avg;

    int maxtime = 0;

    static int autocorr_helper(vector<int> times, int maxtime, int diff);
    static int avg_helper(vector<int> times, int maxtime);
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
        for (long int j = 0; j < Nsites; j++) {
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
    flip_times.resize(system->Nsites);
}

void Measurer::measure()
{
    magnetizations.push_back(system->magnetization);
    energies.push_back(system->energy);
}

void Measurer::spin_flip(int x, int y, int time) {
    flip_times[x*system->Ny + y].push_back(time);
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

int Measurer::autocorr_helper(vector<int> times, int maxtime, int diff)
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

    short flag = +1;

    int curt = diff;
    
    for (int i = 0; (i < times.size()) && (times[i] < diff); i++)
        flag = -flag;
    for (auto p : timesSort) {
        sum += (p - curt)*flag;
        flag = -flag;
        curt = p;
    }

    return sum;
}

int Measurer::avg_helper(vector<int> times, int maxtime)
{
    if (times.size() == 0)
        return maxtime;
    int sum = 0;
    int curt = 0;
    int flag = +1;
    for (int i = 0; i < times.size(); i++) {
        sum += (times[i] - curt)*flag;
        flag = -flag;
        curt = times[i];
    }
    sum += (maxtime - times[times.size()-1])*flag;
    return sum;
}

void Measurer::calc_avg(int maxtime)
{
    spin_avg.resize(system->Nsites);
    for (int i = 0; i < system->Nsites; i++) {
        spin_avg[i] = (double)Measurer::avg_helper(flip_times[i], maxtime) / maxtime;
    }
}

double Measurer::calc_autocorr(int maxtime, int diff)
{
    double sum = 0;
    for (int i = 0; i < system->Nsites; i++) {
        sum += (double)Measurer::autocorr_helper(flip_times[i], maxtime, diff) / (maxtime - diff) - spin_avg[i]*spin_avg[i];
    }
    sum /= system->Nsites;

    return sum;
}


void Measurer::print_results()
{
    cout << "Energy: " << energy << endl;
    cout << "Magnetization: " << magnetization << endl;

}

int main(int argc, char* argv[])
{
    int Nx, Ny, Ttherm, Tmc, corr_time_max;
    double temp, corr_time_dt;
    if (argc != 8) {
        cerr << "Usage: " << argv[0] << " [Nx] [Ny] [temp] [Ttherm] [Tmc] [corr_time_max] [corr_time_dt]" << endl;
        return -1;
    } else {
        Nx = atoi(argv[1]);
        Ny = atoi(argv[2]);
        temp = atof(argv[3]);
        Ttherm = atoi(argv[4]);
        Tmc = atoi(argv[5]);
        corr_time_max = atoi(argv[6]);
        corr_time_dt = atof(argv[7]);
    }

    srand(time(nullptr));

    IsingMCWorker mc(Nx, Ny, 1.0);
    mc.init_lattice_random();

    cerr << "Thermalizing...";
    mc.metropolis(Ttherm, temp);
    cerr << "done" << endl;

    cerr << "Measuring...";

    Measurer m(&mc);
    mc.metropolis(Tmc, temp, &m);

    cerr << "done" << endl;

    cerr << "Calculating...";

    m.calc_results();
    m.calc_avg(Tmc*mc.Nsites);
    m.print_results();

    for (double t = 0; t < corr_time_max; t += corr_time_dt) {
        cout << t << " " << m.calc_autocorr(Tmc*mc.Nsites, round(t*mc.Nsites)) << endl;
    }

    cerr << "done" << endl;

    return 0;
}
