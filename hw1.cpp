#include <stdio.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

using namespace std;

double sum(vector<double> v) {
	double s = 0;
	for (int i = 0; i < v.size(); i++) {
		s += v.at(i);
	}
	return s;
}

double mean(vector<double> v) {
	double m = sum(v) / v.size();
	return m;
}

double median(vector<double> v) {
	vector<double> sor = v;
	sort(sor.begin(), sor.end());
	return sor.at(sor.size()/2);
}

double find_min(vector<double> v) {
	double min = DBL_MAX;
	for (double i : v) {
		if (i < min) {
			min = i;
		}
	}
	return min;
}

double find_max(vector<double> v) {
	double max = DBL_MIN;
	for (double i : v) {
		if (i > max) {
			max = i;
		}
	}
	return max;
}

void range(vector<double> v) {
	cout << find_min(v) << " " << find_max(v);
}

double covar(vector<double> x, vector<double> y) {
	double xmean = mean(x);
	double ymean = mean(y);
	double xsum;
	double ysum;
	double total = 0;
	for (int i = 0; i < x.size(); i++) {
		xsum = x.at(i) - xmean;
		ysum = y.at(i) - ymean;
		total += xsum * ysum;
	}
	return total/(x.size()-1);
}

double cor(vector<double> x, vector<double> y) {
	double xmean = mean(x);
	double ymean = mean(y);
	double xsum = 0;
	double ysum = 0;
	for (int i = 0; i < x.size(); i++) {
		xsum += (x.at(i) - xmean) * (x.at(i) - xmean);
		ysum += (y.at(i) - ymean) * (y.at(i) - ymean);
	}
	double xstd = sqrt(xsum / x.size() );
	double ystd = sqrt(ysum / y.size() );
	return covar(x,y)/(xstd * ystd);
}

int main(int argc, char** argv) {
	ifstream inFS; //Input file stream
	string line;
	string rm_in, medv_in;
	const int MAX_LEN = 1000;
	vector<double> rm(MAX_LEN);
	vector<double> medv(MAX_LEN);

	//Try to open file
	cout << "Opening file Boston.csv." << endl;

	inFS.open("Boston.csv");
	if (!inFS.is_open()) {
		cout << "could not open file Boston.csv" << endl;
		return 1;
	}

	//can now use inFS stream like cin stream
	//Boston.csv should contain two doubles

	cout << "Reading line 1" << endl;
	getline(inFS, line);
	
	//echo heading
	cout << "heading: " << line << endl;

	int numObservations = 0;
	while (inFS.good()) {
		getline(inFS, rm_in, ',');
		getline(inFS, medv_in, '\n');

		rm.at(numObservations) = stof(rm_in);
		medv.at(numObservations) = stof(medv_in);

		numObservations++;
	}

	rm.resize(numObservations);
	medv.resize(numObservations);

	cout << "new length " << rm.size() << endl;

	cout << "Closing file Boston.csv." << endl;
	inFS.close(); //Done wite file, so close it

	cout << "Number of records: " << numObservations << endl;

	cout << "\nStats for rm" << endl;
	//print_stats(rm);
	cout << "Sum: " << sum(rm) << endl;
	cout << "Mean: " << mean(rm) << endl;
	cout << "Median: " << median(rm) << endl;
	cout << "range: ";
	range(rm);
	cout << endl;

	cout << "\nStats for medv" << endl;
	//print_stats(medv);
	cout << "Sum: " << sum(medv) << endl;
	cout << "Mean: " << mean(medv) << endl;
	cout << "Median: " << median(medv) << endl;
	cout << "range: ";
	range(medv);
	cout << endl;

	cout << "\n Covariance = " << covar(rm, medv);
	cout << "\n Correlation = " << cor(rm, medv) << endl;
	cout << "\nProgram terminated.";

	return 0;
}



