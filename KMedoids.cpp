
#include <stdio.h>      
#include <stdlib.h>     
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <math.h>
#include <map>
#include<algorithm>
#include <iomanip> 
#include <float.h> 
#include <limits.h>

using namespace std;
// Class of Packet
class Packet { 
 public:  
	int label;
	string src_addr;
	int src_port;
	string dest_addr;
	int dest_port;
	int protocol;
	int arrival_time;
	int length;
	string addr_port;
};

// Class of Point
class Point { 
 public:  
	string time_str;
	string length_str;
	int medoid;
	int id;
	bool is_medoid;
};

// Set of packets
vector<Packet> packet_vector;
// Set of flows
vector<pair<string, vector<int>>> flows;
// Set of all points
vector<Point> point_vector;
// Set of K medoids
vector<Point> kmedoids;
// Set of points that are not medoids
vector<Point> non_kmedoids;
// Number of initial medoids
int kmedoid_num;
// Set of pairs of selected and non-selected objects
vector<vector<int>> pairs;
// Absolute- error
double absolute_error;

// Read two input files
void ReadFile(char *argv[]) {
	ifstream inf;
	// Read packet information file
    inf.open(argv[1]);  
    if (inf) {
    	 string line;
    	 getline(inf,line);
         while(getline(inf,line)) {
	        stringstream attributes(line);
	        string res;
        	Packet attr;
        	int counter = 0;
            while(attributes>>res) {
            	// Acquire information of the packet
            	if (counter == 0) attr.src_addr = res;
                if (counter == 1) attr.src_port = stoi(res);
                if (counter == 2) attr.dest_addr = res;
                if (counter == 3) attr.dest_port = stoi(res);
                if (counter == 4) attr.protocol = stoi(res);
                if (counter == 5) attr.arrival_time = stoi(res);
                if (counter == 6) attr.length = stoi(res);
                counter++;
            }
            attr.addr_port = attr.src_addr+to_string(attr.src_port)+attr.dest_addr+ to_string(attr.dest_port);
			packet_vector.push_back(attr);
        }
    }
    inf.close();
    
    // Read initial set of K medoid file
	inf.open(argv[2]); 
	if (inf) {
		string line;
		getline(inf,line);
    	kmedoid_num = stoi(line);
    	getline(inf,line);
    	stringstream attributes(line);
    	Point pt;
    	string res;
    	for (int i = 0; i < kmedoid_num; ++i) {
    		attributes>>res;
    		// Acquire information of the K medoid
    		pt.id = stoi(res);
    		pt.medoid = stoi(res);
    		kmedoids.push_back(pt);
    	}
    }
    inf.close();
}

void PrintFlow() {
    ofstream out1( "Flow.txt" );
    int id = 0;
    for (int i = 0; i < flows.size(); ++i) {
    	if (((flows[i].second).size() > 1)) {
    		// Calculate the average transferring time and the average packet length of the flow.
    		double avg_time = 0;
    		double avg_length = 0;
    		vector<int> cur_indices = flows[i].second;
    		for (int j = cur_indices.size() - 1; j >= 0; --j) {
    			if (j > 0) {
    				avg_time += packet_vector[cur_indices[j]].arrival_time - packet_vector[cur_indices[j-1]].arrival_time;
    			}
    			avg_length += packet_vector[cur_indices[j]].length;
    		}
    		avg_time = avg_time / (cur_indices.size() - 1);
    		avg_length = avg_length / cur_indices.size();

    		// Print flow data to Flow.txt
    		out1 << id << " " << fixed << setprecision(2) << avg_time << " " << fixed << setprecision(2) << avg_length<<endl;
    		
    		// Store information of the flow
    		Point pt;
    		pt.time_str = to_string(avg_time);
    		pt.length_str = to_string(avg_length);
    		pt.id = id;
    		pt.is_medoid = false;
    		point_vector.push_back(pt);
    		id++;
    	}
    }  
    out1.close();	
};

// Data reprocessing module
void DataReprocessing(char *argv[]) {
	// Read two input files
	ReadFile(argv);
 	// Merge the packets into flows
    for (int i = 0; i < packet_vector.size(); ++i){
    	vector<int> indices;   
    	int flag = -1;
    	for (int j = 0; j < flows.size(); ++j) {
    		if (flows[j].first == packet_vector[i].addr_port) {
    			flag = j;
    		}
    	}
    	if(flag != -1) {
    		vector<int> cur_indices = flows[flag].second;
    		cur_indices.push_back(i);
			flows[flag].second = cur_indices;
    	} else {
    		indices.push_back(i);
    		pair<string, vector<int>> p1(packet_vector[i].addr_port, indices);

    		flows.push_back(p1);
    	}
    	indices.clear();
    }
    PrintFlow();
};

// Sort flows by ID
bool SortById(Point a, Point b) { 
    return (a.id < b.id); 
};

// Calculate absolute-error
double CalculateAbsoluteError() {
	double res = 0;
	for (int i = 0; i < non_kmedoids.size(); ++i) {
		int min_pt = -1;
		int min_pt_in_med_vec = -1;
		double min_dissimilarity = 10000000000;
		for (int j = 0; j < kmedoid_num; ++j) {
			// Use Mannhaton distance to measure the distance between flow
			double cal_res =  abs( stod(non_kmedoids[i].time_str) - stod(kmedoids[j].time_str) )  + abs(stod(non_kmedoids[i].length_str) - stod(kmedoids[j].length_str));
			double tcost = cal_res - min_dissimilarity;
			if (tcost < 0 && fabs(tcost)> 0.001){
				min_dissimilarity = cal_res;
				min_pt = kmedoids[j].id;
				min_pt_in_med_vec = j;
			}
		}
		non_kmedoids[i].medoid = min_pt;
		res+=min_dissimilarity;
	}
	return res;
};

// Calculate total swapping cost
bool CalculateTotalSwappingCost(bool flag) {
		// Record current K medoids
	vector<int> med_indices;
	for (int a = 0; a < kmedoids.size(); ++a) {
		med_indices.push_back(kmedoids[a].id);
	}
	
	
	bool break_loops = false;
	for (int i = 0; i < non_kmedoids.size(); ++i) {
		
		if (break_loops) break;
		
		Point this_non_medoid = non_kmedoids[i];
		for (int j = 0; j < kmedoids.size(); ++j) {
			// Use a flow as a new mediod
			Point this_medoid = kmedoids[j];

			bool found_duplicate = false;
			int asdc = 0;
			for (int d = 0; d < pairs.size(); ++d) {
				vector<int> current_pair = pairs[d];
				if ((current_pair[0] == this_medoid.id && current_pair[1] == this_non_medoid.id) || (current_pair[1] == this_medoid.id && current_pair[0] == this_non_medoid.id)) {
					found_duplicate = true;
				}
			}
			// If the pair of selected and non-selected objects is already in the set, skips these pair
			if (found_duplicate) continue;

			kmedoids[j] = this_non_medoid;
			non_kmedoids[i] = this_medoid;
			double currentAbsoluteError = CalculateAbsoluteError();
			double total_swapping_cost = currentAbsoluteError - absolute_error;

			// If the total swapping cost is smaller than zero, the swaps is carried out
			if (total_swapping_cost < 0 && abs(total_swapping_cost)> 0.001){
				absolute_error = currentAbsoluteError;
				break_loops = true;
				vector<int> points;
				points.push_back(this_medoid.id);
				points.push_back(this_non_medoid.id);
				pairs.push_back(points);
				break;
			} else {
				kmedoids[j] = this_medoid;
				non_kmedoids[i] = this_non_medoid;
			}
		}
	}
	int count_med = 0;
	for (int b = 0; b < kmedoids.size(); ++b) {
		for (int c = 0; c < med_indices.size(); ++c) {
			if (med_indices[c] == kmedoids[b].id) {
				++count_med;
			}
		}
	}
	// If there is no change of the medoids, ends the while loop
	if (count_med == kmedoids.size()) flag = true;
	return flag;
}

void PrintKMedoidsClusters() {
	absolute_error = CalculateAbsoluteError();
	ofstream out2( "KMedoidsClusters.txt" );
	out2 << fixed << setprecision(2) << absolute_error<<endl;
	for (int i = 0; i < kmedoids.size(); ++i) {
		out2<<kmedoids[i].id;
		if (i != kmedoids.size()-1) {
			out2<<" ";
		}
	}
	out2<<endl;
	for (int i = 0; i < kmedoids.size(); ++i) {
		vector<Point> id_vec;
		for (int j = 0; j < non_kmedoids.size(); ++j) {
			if (non_kmedoids[j].medoid == kmedoids[i].id) {
				id_vec.push_back(non_kmedoids[j]);
			}
		}
		id_vec.push_back(kmedoids[i]);
		std::sort(id_vec.begin(), id_vec.end(), SortById);
		for (int k = 0; k < id_vec.size(); ++k) {
			out2<<id_vec[k].id;
			if (k != id_vec.size()-1) {
				out2<<" ";
			}
		} 
		out2<<endl;
	}
    out2.close();
};

// Clustering module
void Clustering() {
	// Prepare the set of medoids
	for (int i = 0; i < point_vector.size(); ++i) {
		for (int j = 0; j < kmedoid_num; ++j) {
			if (i == kmedoids[j].id ) {
				 point_vector[i].medoid = i;	
				 point_vector[i].is_medoid = true;
				 kmedoids[j].time_str = point_vector[i].time_str;
				 kmedoids[j].length_str = point_vector[i].length_str;
			}
		}
	}
	 // Prepare the set of points that are not medoids
	for (int i = 0; i < point_vector.size(); ++i) {
		if (!point_vector[i].is_medoid) {
			non_kmedoids.push_back(point_vector[i]);
		}
	}

	bool stop_flag = false;
	while(!stop_flag) {
		absolute_error = CalculateAbsoluteError();
		stop_flag = CalculateTotalSwappingCost(stop_flag);
	}
	for (int i = 0; i < kmedoids.size(); ++i) {
		vector<Point> id_vec;
		int second_pos = 0;
		for (int j = 0; j < non_kmedoids.size(); ++j) {
			if (non_kmedoids[j].medoid == kmedoids[i].id) {
				id_vec.push_back(non_kmedoids[j]);
				second_pos = j;
			}
		}
		id_vec.push_back(kmedoids[i]);
		std::sort(id_vec.begin(), id_vec.end(), SortById);
		// If there are two points in a cluster, uses the one with smaller as the medoid
		if (id_vec.size() == 2) {
			kmedoids[i] = id_vec[0];
			non_kmedoids[second_pos] = id_vec[1];
		}
	}
	PrintKMedoidsClusters();
}

int main(int argc,char *argv[])
{
    DataReprocessing(argv);
    Clustering();
    return 0;
}
