#include <iostream>
#include <random>
#include <locale>
#include <fstream>
#include <map>
#include <string> 

using namespace std;

#define SIZE 100

short lattice[SIZE][SIZE];
int steps;
float T;
double w[5];
double M, E;

vector <double> corr, singleCorr, singleSpinCorr, Ham, Mom;
vector <int> singleMom[SIZE][SIZE];
vector <short> singleSpin;

int n = SIZE * SIZE;


void calcW()
{
	// exp array

	double e4 = exp(-4 / T), e8 = e4 * e4;
	w[0] = w[4] = e8;
	w[1] = w[3] = e4;
	w[2] = 0;
}

void set_lattice()
{
	// sets lattice of random chosen spins

	M = E = 0;
	static std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0, RAND_MAX);

	for (int i = 0; i < SIZE; i++)
		for (int j = 0; j < SIZE; j++)
		{
			lattice[i][j] =  ((distribution(generator) / (double)RAND_MAX >= 0.5) - 1) * 2 + 1;
			M += lattice[i][j];
		}

	for (int i = 0; i < SIZE; i++)
		for (int j = 0; j < SIZE; j++)
		{
			E += (i + 1 != SIZE) ? lattice[i][j] * lattice[i + 1][j] : 0;
			E += (j + 1 != SIZE) ? lattice[i][j] * lattice[i][j + 1] : 0;
		}

	calcW();
}

void metropolis()
{
	// update lattice, energy, magnetization

	int x, y, sum;
	for (int i = 0; i < steps; i++)
	{
		x = rand() % SIZE;
		y = rand() % SIZE;
		sum = lattice[(x - 1 + SIZE) % SIZE][y] +
			lattice[(x + 1 + SIZE) % SIZE][y] +
			lattice[x][(y - 1 + SIZE) % SIZE] +
			lattice[x][(y + 1 + SIZE) % SIZE];


		if (sum * lattice[x][y] <= 0 || (rand() / (double)RAND_MAX) < w[sum / 2 + 2])
		{
			lattice[x][y] = -lattice[x][y];
			M += 2 * lattice[x][y];
			E -= 2 * lattice[x][y] * sum;

			singleMom[x][y].push_back(i);
		}

		Ham.push_back(E);
		Mom.push_back(M);
	}
}


void autocorr()
{
	// magnetization autocorrelation function 

	double sum = 0;

	for (int i = 0; i < steps; i ++) {
		sum = 0;
		for (int j = i; j < steps; j++) {
			sum += Mom[j] * Mom[j - i];
		}
		corr.push_back(sum / (steps - i));
		i += 9999;
	}
}

void recover_spin(vector <int> times)
{
	// recovers single spin evolution from spin flip moments

	int len = times.size();
	short s = 1;

	singleSpin.clear();

	if (len == 0) {
		for (int count = 0; count < steps; count++) {
			singleSpin.push_back(s);
		}
	}
	else
	{
		for (int count = 0; count < times[0]; count++) {
			singleSpin.push_back(s);
		}

		for (int k = 1; k < len; k++) {
			s = -s;
			for (int count = 0; count < (times[k] - times[k - 1]); count++) {
				singleSpin.push_back(s);
			}
		}

		s = -s;

		for (int count = 0; count < (steps - times[len - 1]); count++) {
			singleSpin.push_back(s);
		}

	}
}

int prod(vector <short> vec, int diff)
{
	// returns scalar product of a given vector with its index-shifted copy

	int len = vec.size();
	int sum = 0;

	for (int count = 0; count < len - diff; count++) {
		sum += vec[count] * vec[count + diff];
	}

	return sum;
}

vector <double> mean_vec(vector <double> vec1, int weight1, vector <double> vec2, int weight2) 
{
	//  returns weighted mean vector

	vector <double> meanVec;
	int len = vec1.size();
	int weight = weight1 + weight2;

	for (int count = 0; count < len; count++) {
		meanVec.push_back((vec1[count] * weight1 + vec2[count] * weight2) / weight);
	}

	return meanVec;
}

int dict_sum_d(vector <int> times, vector <int> timesd, int diff)
{
	// counts autocorr

	int sum = 0, len = times.size(), prev = diff, lenSort, min = 0;
	short flag = 1;
	map<int, int> dict;
	vector<int> timesSort;
	
	for (int i = 0; i < len; i++) {
		if (times[i] > diff)
			dict[times[i]] += (i + 1);
		if (timesd[i] < steps - diff)
			dict[timesd[i]] -= (i - 1);
	}
	
	for (int i = 0; i < len; i++) {
		if (times[i] > diff)
			timesSort.push_back(times[i]);
		if (timesd[i] < steps - diff)
			timesSort.push_back(timesd[i]);
	}
			
	timesSort.push_back(diff);
	timesSort.push_back(steps - diff);
	sort(timesSort.begin(), timesSort.end());
	lenSort = timesSort.size();

	for (int i = 0; i < len; i++)
		if (times[i] - diff > 0)
			if (dict[times[i]] % 2 != 0)
				flag = -1;

	for (int count = 1; count < lenSort; count++) {
		if (dict[timesSort[count]] != 0) {
			sum += flag * (timesSort[count] - timesSort[count - 1]);
			flag = -flag;
		}
		else {
			sum += flag * (timesSort[count] - timesSort[count - 1]);
		}
	}

	return sum;
}

vector <int> shift_vec_values(vector <int> times, int diff)
{
	// adds constant to vector values

	int len = times.size();
	vector <int> timesd;

	for (int i = 0; i < len; i++) {
		timesd.push_back(times[i] + diff);
	}

	return timesd;
}

void autocorr_single()
{
	// mean single spin autocorrelation S_ii

	int spinNum = 0;
	
	// autocorrelation func of s_00
	for (int d = 0; d < steps; d ++) {
		singleCorr.push_back(dict_sum_d(singleMom[0][0], shift_vec_values(singleMom[0][0], d), d));
		d += 9999;
	}

	// averaging autocorrelation funcs of s_ji

	for (int i = 0; i < SIZE; i++) {
		for (int j = 0; j < SIZE; j++) {
			if (i == 0 && j == 0) continue;
			spinNum++;
			singleSpinCorr.clear();
			for (int d = 0; d < steps; d++) {
				singleSpinCorr.push_back(dict_sum_d(singleMom[i][j], shift_vec_values(singleMom[i][j], d), d));
				d += 9999;
			}
			singleCorr = mean_vec(singleCorr, spinNum, singleSpinCorr, 1);
		}
	}
}

int main()
{
	cout << "Temp:";
	cin >> T;
	cout << "num of spin flips:";
	cin >> steps;

	set_lattice();
	cout << endl << "2d square lattice with size " << SIZE << " is set" << endl;

	metropolis();

	cout << "starting autocorr" << endl;
	autocorr();
	cout << "autocorr done" << endl;

	cout << "starting single autocorr" << endl;
	autocorr_single();
	cout << "single autocorr done" << endl;

	ofstream myfile("ham" + std::to_string(T) + ".txt");
	if (myfile.is_open())
	{
		for (int count = 0; count < steps; count++) {
			myfile << Ham[count] << " ";
		}
		myfile.close();
	}
	else cout << "Unable to open file";

	ofstream myfile1("mom" + std::to_string(T) + ".txt");
	if (myfile1.is_open())
	{
		for (int count1 = 0; count1 < steps; count1++) {
			myfile1 << Mom[count1] << " ";
		}
		myfile1.close();
	}
	else cout << "Unable to open file";

	ofstream myfile2("corr" + std::to_string(T) + ".txt");
	if (myfile2.is_open())
	{
		for (int count1 = 0; count1 < corr.size(); count1++) {
			myfile2 << corr[count1] << " ";
		}
		myfile2.close();
	}
	else cout << "Unable to open file";


	ofstream myfile3("single_corr" + std::to_string(T) + ".txt");
	if (myfile3.is_open())
	{
		for (int count1 = 0; count1 < singleCorr.size(); count1++) {
			myfile3 << singleCorr[count1] << " ";
		}
		myfile3.close();
	}
	else cout << "Unable to open file";

}
