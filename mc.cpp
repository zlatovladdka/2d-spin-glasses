#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <map>
#include <numeric>
#include <random>
#include <stdio.h>

#include "argh.h"

using namespace std;

static random_device r;
static default_random_engine engine(r());

class Measurer;

class IsingMCWorker {
public:
    IsingMCWorker(int _Nx, int _Ny, double _J);

    void init_J_lattice_random(long long seed, int distType);
    void init_J_lattice_const();
    void init_lattice_random();
    void print_J_lattice();
    void print_system();

    inline double& Jlatv(int x, int y) {
        // Periodic boundary conditions
        return jLatticeV[((x + Nx) % Nx) * Ny + (y + Ny) % Ny];
    }

    inline double& Jlath(int x, int y) {
        // Periodic boundary conditions
        return jLatticeH[((x + Nx) % Nx) * Ny + (y + Ny) % Ny];
    }

    inline short& lat(int x, int y) {
        // Periodic boundary conditions
        return lattice[((x + Nx)%Nx)*Ny + (y+Ny)%Ny];
    }

    void metropolis(long long Nsteps, double T, Measurer *m);

    // Lattice parameters
    int Nx, Ny, Nsites;
    vector<short> lattice;
    vector<double> jLatticeV, jLatticeH;
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

    void print_spin_average();
    void print_energy_history();
    void print_magnetization_history();

    void spin_flip(int x, int y, long long time);
    void calc_avg(long long maxtime);
    double calc_autocorr(long long maxtime, long long diff);

    double energy = 0.0;
    double energy_fluct = 0.0;
    double magnetization = 0.0;
    double magnetization_fluct = 0.0;
private:
    IsingMCWorker* system;

    vector<int> magnetizations;
    vector<double> energies;

    vector<vector<long long>> flip_times;
    vector<double> spin_avg;

    long long maxtime = 0;

    static long long autocorr_helper(vector<long long> times, long long maxtime, long long diff);
    static long long avg_helper(vector<long long> times, long long maxtime);
};

IsingMCWorker::IsingMCWorker(int _Nx, int _Ny, double _J)
    : Nx(_Nx), Ny(_Ny), Nsites(Nx*Ny), J(_J)
{
    energy = 0.0;
    magnetization = 0.0;
    lattice.resize(Nsites);
    jLatticeV.resize(Nsites);
    jLatticeH.resize(Nsites);
}

void IsingMCWorker::init_lattice_random()
{
    uniform_int_distribution<int> spin_dist(0, 1);
    magnetization = 0;
    generate(lattice.begin(), lattice.end(), [&]() {
                short s = 2 * spin_dist(engine) - 1;
                magnetization += s;
                return s;
            });

    energy = 0.0;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            // energy += -J*lat(i, j)*(lat(i, j+1) + lat(i+1,j));

            energy += -J * lat(i, j) * (Jlath(i, j) * lat(i, j + 1) + Jlatv(i, j) * lat(i + 1, j));
        }
    }
}

void IsingMCWorker::init_J_lattice_random(long long seed, int distType)
{
    static random_device jr;
    static default_random_engine jengine(jr());
    jengine.seed(seed);

    if (distType == 0) {
        normal_distribution<double> jStrength(0, J);

        generate(jLatticeV.begin(), jLatticeV.end(), [&]() {
            return jStrength(jengine);
            });

        generate(jLatticeH.begin(), jLatticeH.end(), [&]() {
            return jStrength(jengine);
            });
    }
    else {
        uniform_int_distribution<> jStrength(0, 1);

        generate(jLatticeV.begin(), jLatticeV.end(), [&]() {
            return J * (jStrength(jengine) * 2 - 1);
            });

        generate(jLatticeH.begin(), jLatticeH.end(), [&]() {
            return J * (jStrength(jengine) * 2 - 1);
            });
    }
}

void IsingMCWorker::init_J_lattice_const()
{
    fill(jLatticeV.begin(), jLatticeV.end(), J);
    fill(jLatticeH.begin(), jLatticeH.end(), J);
}

void IsingMCWorker::print_system() {
    cerr << " === System state === " << endl;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            cerr << (lat(i,j) == +1 ? '+' : '-');
        }
        cerr << endl;
    }
    cerr << "Magnetization: " << magnetization << endl;
    cerr << "Energy: " << energy << endl;
}

void IsingMCWorker::print_J_lattice() {
    cerr << " === J-lattice (vertical) === " << endl;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            cerr << Jlatv(i, j) << " ";
        }
        cerr << endl;
    }

    cerr << " === J-lattice (horizontal)=== " << endl;
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            cerr << Jlath(i, j) << " ";
        }
        cerr << endl;
    }
}

void IsingMCWorker::metropolis(long long time, double T, Measurer *m = nullptr)
{
    uniform_int_distribution<int> distx(0, Nx-1);
    uniform_int_distribution<int> disty(0, Ny-1);
    uniform_real_distribution<double> distmc(0, 1);
	for (long long i = 0; i < time; i++)
	{
        for (int j = 0; j < Nsites; j++) {
            int x = distx(engine);
            int y = disty(engine);
            // int sum = lat(x-1,y) + lat(x+1, y) + lat(x,y-1) + lat(x, y+1);

            double sum = lat(x-1,y)*Jlatv(x-1,y)+lat(x+1,y)*Jlatv(x, y)+lat(x, y - 1)*Jlath(x,y-1)+lat(x, y + 1)*Jlath(x,y);
            double dE = 2*lat(x,y)*sum*J;

            if (dE <= 0 || distmc(engine) < exp(-dE/T))
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

void Measurer::spin_flip(int x, int y, long long time) {
    flip_times[x*system->Ny + y].push_back(time);
}

void Measurer::calc_results()
{
    energy = accumulate(energies.begin(), energies.end(), 0.0);
    energy /= energies.size();
    energy /= system->Nsites;

    energy_fluct = accumulate(energies.begin(), energies.end(), 0.0,
            [](double sum, double energy){ return sum + energy*energy; } );
    energy_fluct /= energies.size();
    energy_fluct /= system->Nsites*system->Nsites;
    energy_fluct -= energy*energy;
    energy_fluct *= system->Nsites;

    magnetization = accumulate(magnetizations.begin(), magnetizations.end(), 0.0);
    magnetization /= magnetizations.size();
    magnetization /= system->Nsites;

    magnetization_fluct = accumulate(magnetizations.begin(), magnetizations.end(), 0.0,
            [](double sum, double magn){ return sum + magn*magn; } );
    magnetization_fluct /= magnetizations.size();
    magnetization_fluct /= system->Nsites*system->Nsites;
    magnetization_fluct -= magnetization*magnetization;
    magnetization_fluct *= system->Nsites;

}

long long Measurer::autocorr_helper(vector<long long> times, long long maxtime, long long diff)
{
    long long sum = 0;
    vector<long long> timesSort;

    for (long long i = 0; i < times.size(); i++) {
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

    long long curt = diff;

    for (long long i = 0; (i < times.size()) && (times[i] < diff); i++)
        flag = -flag;
    for (auto p : timesSort) {
        sum += (p - curt)*flag;
        flag = -flag;
        curt = p;
    }

    return sum;
}

long long Measurer::avg_helper(vector<long long> times, long long maxtime)
{
    if (times.size() == 0)
        return maxtime;
    long long sum = 0;
    long long curt = 0;
    int flag = +1;
    for (long long i = 0; i < times.size(); i++) {
        sum += (times[i] - curt)*flag;
        flag = -flag;
        curt = times[i];
    }
    sum += (maxtime - curt)*flag;
    return sum;
}

void Measurer::calc_avg(long long maxtime)
{
    spin_avg.resize(system->Nsites);
    for (int i = 0; i < system->Nsites; i++) {
        spin_avg[i] = (double)Measurer::avg_helper(flip_times[i], maxtime) / maxtime;
    }
}

double Measurer::calc_autocorr(long long maxtime, long long diff)
{
    double sum = 0;
    for (int i = 0; i < system->Nsites; i++) {
        sum += (double)Measurer::autocorr_helper(flip_times[i], maxtime, diff) / (maxtime - diff) - spin_avg[i]*spin_avg[i];
    }
    sum /= system->Nsites;

    return sum;
}

void Measurer::print_spin_average() {
    for (int i = 0; i < system->Nx; i++) {
        for (int j = 0; j < system->Ny; j++) {
            cout << spin_avg[i*system->Ny + j] << " ";
        }
        cout << endl;
    }
}

void Measurer::print_energy_history() {
    for (int i = 0; i < energies.size(); i++) {
        cout << " " << energies[i];
    }
}

void Measurer::print_magnetization_history() {
    for (int i = 0; i < magnetizations.size(); i++) {
        cout << i << " " << magnetizations[i] << endl;
    }
}

void Measurer::print_results()
{
    cerr << "Energy: ";
    cout << energy << " " << energy_fluct << endl;
    cerr << "Magnetization: ";
    cout << magnetization << " " << magnetization_fluct << endl;

}

static void print_usage(argh::parser& cmd) {
    cerr << "Usage: " << cmd(0).str() << " [OPTIONS]" << endl;
    cerr << endl;
    cerr << "Available options are:" << endl;
    cerr << "   -v, --verbose   Be more verbose" << endl;
    cerr << "   -h, --help      Print this help message" << endl;
    cerr << endl;
    cerr << "System parameters:" << endl;
    cerr << "   -x, --nx NX     System size along X direction" << endl;
    cerr << "   -y, --ny NY     System size along Y direction" << endl;
    cerr << "   --temp          Temperature for MC simulation" << endl;
    cerr << "   -J J            Coupling constant (default: 1.0)" << endl;
    cerr << "   -g, --glass     Initialize random J couplings (default: constant)" << endl;
    cerr << "   -d, --dist      Distribution (0 - Gaussian, 1 - plus-minus J) (default: 0 - Gaussian)" << endl;
    cerr << "   --therm TIME    Run TIME thermalization sweeps before measuring (default: 0)" << endl;
    cerr << "   --time TIME     Run TIME MC sweeps with measurements" << endl;
    cerr << "   --seed SEED     Random seed for J-lattice generation" << endl;
    cerr << endl;
    cerr << "Results reporting:" << endl;
    cerr << "   --spin-average              Print average values of each spin" << endl;
    cerr << "   --energy-history            Print energy dependence on the MC time" << endl;
    cerr << "   --magnetization-history     Print magnetization dependence on the MC time" << endl;
    cerr << "   --autocorr MAXTIME          Calculate single-spin autocorrelation function up to time MAXTIME" << endl;
    cerr << "   --autocorr-dt DT            Calculate autocorrelation function each DT sweeps (default: 1.0)" << endl;
}

int main(int argc, char* argv[])
{
    // 
    // Command line parsing
    //
    argh::parser cmd(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);

    if (cmd[{"-h", "--help"}]) {
        print_usage(cmd);
        return 0;
    }
    bool verbose = cmd[ {"-v", "--verbose"} ];
    // System size
    int Nx, Ny;
    if (!(cmd({"-x", "--nx"}) >> Nx) || !(cmd({"-y", "--ny"}) >> Ny)) {
        cerr << "Please specify system size!" << endl;
        print_usage(cmd);
        return 1;
    }
    // Coupling constant
    double J;
    cmd("-J", 1.0) >> J;

    // Coupling distribution
    int dist;
    if (!(cmd({ "-d", "--dist" }) >> dist)) {
        dist = 0;
    }

    // Initializing MC worker
    IsingMCWorker mc(Nx, Ny, J);

    long long randSeed;
    if (cmd[{"-g", "--glass"}]) {
        if (!(cmd({ "--seed" }) >> randSeed)) {
            cerr << "Please specify random seed!" << endl;
            print_usage(cmd);
            return 1;
        }
        cerr << "Initializing random couplings" << endl;
        mc.init_J_lattice_random(randSeed, dist);   
    } else {
        cerr << "Initializing constant couplings" << endl;
        mc.init_J_lattice_const();
    }
    mc.init_lattice_random();

    if (verbose) {
        mc.print_J_lattice();
    }

    // Temperature
    double temp;
    if (!(cmd({"-t", "--temp"}) >> temp)) {
        cerr << "Please specify temperature!" << endl;
        print_usage(cmd);
        return 1;
    }

    // Thermalization
    long long Ttherm;
    cmd("--therm", 0) >> Ttherm;
    if (Ttherm > 0) {
        cerr << "Thermalizing... ";
        mc.metropolis(Ttherm, temp);
        cerr << "done!" << endl;
    }

    if (verbose) {
        mc.print_system();
    }
    // Measuring
    long long Tmc;
    if (!(cmd("--time") >> Tmc)) {
        cerr << "Please provide Monte-Carlo time!" << endl;
        print_usage(cmd);
        return 1;
    }
    cerr << "Running Monte-Carlo... ";
    Measurer m(&mc);
    mc.metropolis(Tmc, temp, &m);
    cerr << "done!" << endl;

    m.calc_results();
    m.calc_avg(Tmc*mc.Nsites);

    m.print_results();

    if (cmd["--spin-average"]) {
        cerr << "Average values of spins:" << endl;
        m.print_spin_average();
    }

    if (cmd["--energy-history"]) {
        cerr << "Energies:" << endl;
        m.print_energy_history();
    }
    if (cmd["--magnetization-history"]) {
        cerr << "Magnetizations:" << endl;
        m.print_magnetization_history();
    }
    double autocorr, autocorr_dt;
    if ((cmd("--autocorr") >> autocorr) && (cmd("--autocorr-dt", 1.0) >> autocorr_dt)) {
        cerr << "Single-spin autocorrelations:" << endl;
        for (double t = 0; t < autocorr; t += autocorr_dt) {
            cout << t << " " << m.calc_autocorr(Tmc*mc.Nsites, round(t*mc.Nsites)) << endl;
        }
    }

    return 0;
}
