# TrafficPacketClusteringEngine
A traffic packet clustering engine to cluster the raw network packet to different applications, such as http, smtp.

Data preprocessing module and clustering module have been implemented.

Using following command:
g++ -std=c++11 KMedoids.cpp -o runKMedoids
./runKMedoids network packets.txt initial medoids.txt
