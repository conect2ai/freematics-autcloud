#include <iostream>
#include <vector>
#include <set>
#include <cmath>

// define global file to save the results
// std::ofstream outfile ("debug.csv");
// std::ofstream outfile ("output.csv");

class DataCloud {
public:
    // Constructor
    DataCloud(const std::vector<double>& x, int point_key) {
        n = 1;
        points.push_back(point_key);
        mean = x;
        var = 0;
        pertinency = 1;
    }

    // External methods
    std::vector<double> getUpdatedMean(const std::vector<double>& x) {
        int new_n = n + 1;
        std::vector<double> updatedMean(mean.size());
        for (size_t i = 0; i < mean.size(); ++i) {
            updatedMean[i] = ((new_n - 1) * mean[i] + x[i]) / new_n;
        }
        return updatedMean;
    }

    double getUpdatedVarianceK2(const std::vector<double>& x) {
        int new_n = n + 1;
        double distanceSquared = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            double diff = x[i] - mean[i];
            distanceSquared += diff * diff;
        }
        return std::sqrt(distanceSquared) / 2.0;
    }

    double getUpdatedVarianceK3(const std::vector<double>& x, const std::vector<double>& test_mean) {
        int new_n = n + 1;
        double distanceSquared = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            double diff = x[i] - test_mean[i];
            distanceSquared += diff * diff;
        }
        return ((new_n - 1) * var + distanceSquared) / new_n;
    }

    void updateStats(const std::vector<double>& new_mean, double new_var) {
        mean = new_mean;
        var = new_var;
    }

    void addPoint(int new_point) {
        points.push_back(new_point);
        n += 1;
    }

    void addPoints(const std::vector<int>& new_points) {
        points.insert(points.end(), new_points.begin(), new_points.end());
        // n = points.size();
        // remove duplicates
        std::sort(points.begin(), points.end());
        points.erase(std::unique(points.begin(), points.end()), points.end());
        // n += new_points.size();
        n = points.size();
        
    }

    double calculate_force(const double point) const {
        // Calcula a 'força' baseada na distância do ponto à média
        double distance = std::abs(point - mean[0]);  // Assumindo que mean é um vetor de uma dimensão
        // Em ForceAtlas, a força geralmente diminui com a distância
        // Você pode experimentar com diferentes funções aqui
        double force = 1.0 / (1.0 + distance);
        return force;
    }

    void adjust_variance_with_force(const std::vector<double>& data) {
        // Ajustar a variância com base na força
        double total_force = 0.0;
        for (const auto& point : data) {
            total_force += calculate_force(point);
        }

        // Ajustar a variância com base na força total
        // Substitua com a implementação real
        var += total_force;
    }

    std::vector<int> getPoints() const {
        return points;
    }

    int getN() {
        return n;
    }

    bool operator==(const DataCloud& c) {
        return points == c.points;
    }

public:
    int n;
    std::vector<int> points;
    std::vector<double> mean;
    double var;
    double pertinency;
};

class TEDACloud {
public:
    TEDACloud() {
        index_point = 0;
        classIndex = {{1.0}, {1.0}};
        alfa = {0.0, 0.0};
        relevanceList = {0.0, 0.0};
    }

    void createCloud(const std::vector<double>& x, int k) {
        DataCloud c(x, k);
        clouds.push_back(c);
    }

    void deleteCloud(DataCloud& c) {
        for (size_t i = 0; i < clouds.size(); ++i) {
            if (clouds[i] == c) {
                clouds.erase(clouds.begin() + i);
                break;
            }
        }
    }

    DataCloud& getCloud(int index) {
        return clouds[index];
    }

    int getNumClouds() {
        // return clouds.size();
        return 3;
    }


    void print_cloud_points() {
        int totalPoints = 0;
        // iterate over clouds and print the quantity of points in each cloud
        for (size_t i = 0; i < clouds.size(); ++i) {
            std::cout << "Cloud " << i << " has " << clouds[i].getPoints().size() << " points." << std::endl;
            totalPoints += clouds[i].getPoints().size();
        }
        std::cout << "Total points: " << totalPoints << std::endl;
    }

    void print_cloud_index_points() {
        // print the points index of each cloud
        for (size_t i = 0; i < clouds.size(); ++i) {
            std::cout << "Cloud " << i << " has the following points: ";
            for (size_t j = 0; j < clouds[i].getPoints().size(); ++j) {
                std::cout << clouds[i].getPoints()[j] << ", ";
            }
            std::cout << std::endl;
        }
    }

    void metrics() {
        std::cout << "Number of DataClouds: " << clouds.size() << std::endl;
    }

    int runOnline(const std::vector<double>& data, int m, int numCloud = 3, bool isOutlier = false) {
        
        index_point++; // Atualize o contador de linhas

        if (index_point == 1) {
            // Crie a primeira nuvem de dados e adicione X
            createCloud(data, index_point);
            // printf("Criando a primeira nuvem de dados.\n");
            return 0;
        }
        if (index_point == 2){
            createCloud(data, index_point);
            // printf("Criando a segunda nuvem de dados.\n");
            return 1;
        } 
        if (index_point >= 3){
            
            // Verifique se o vetor de clouds está cheio
            if (index_point == 3){
                createCloud(data, index_point);
                // std::cout << "Criando nova nuvem para o ponto: " << index_point << " - Quantidade de nuvens: " << clouds.size() + 1 << std::endl;
                return 2;
            }

            bool shouldCreateNewCloud = true;
            float maxTipycality = -std::numeric_limits<float>::infinity();
            int chosenCloud = -1;
            bool nothingCloud = true;
            // alfa = std::vector<double>(getNumClouds(), 0.0);

            for (int i = 0; i < numCloud; i++) {
                DataCloud& c = getCloud(i);

                // Calcule as novas médias e variâncias se adicionássemos o ponto
                int test_n = c.getN() + 1;
                std::vector<double> test_mean = c.getUpdatedMean(data);
                double test_var = c.getUpdatedVarianceK3(data, test_mean);

                std::vector<double> diff(test_mean.size());
                double sum_diff_squared = 0.0;

                for (size_t j = 0; j < test_mean.size(); ++j) {
                    diff[j] = test_mean[j] - data[j];
                    double diff_squared = diff[j] * diff[j];
                    sum_diff_squared += diff_squared;
                }

                double eccentricity = (test_var + sum_diff_squared) / (test_n * test_var);

                // Escreva a excentricidade no arquivo
                // outfile << "Point: " << index_point << " eccentricity: " << eccentricity << std::endl;

                double norm_eccentricity = eccentricity / 2.0;
                // Calcule a tipicidade e tipicidade normalizada
                double typicality = 1 - eccentricity;
                double norm_typicality = typicality / (index_point - 2);

                // outfile << "Point: " << index_point << " norm_typicality: " << norm_typicality << std::endl;

                bool is_point_in_c = norm_eccentricity <= (std::pow(m, 2) + 1) / (2 * test_n);

                if (is_point_in_c) {
                    // Atualize as estatísticas da nuvem
                    c.updateStats(test_mean, test_var);
                    // Adicione o ponto à nuvem
                    c.addPoint(index_point);
                    shouldCreateNewCloud = false;
                    nothingCloud = false;
                    chosenCloud = i;
                } else if (norm_typicality >= maxTipycality && nothingCloud) {
                    maxTipycality = norm_typicality;
                    chosenCloud = i;
                    shouldCreateNewCloud = false;
                }

                alfa[i] = norm_typicality;
            }

            // if (nothingCloud && chosenCloud != -1){
                
            //     if (isOutlier){
            //         // get chosen cloud
            //         DataCloud& c = getCloud(chosenCloud);

            //         c.mean = data;

            //         c.adjust_variance_with_force(data);

            //         Serial.println("Clouds reinitialized");

            //     } else {
            //         std::vector<double> test_mean = getCloud(chosenCloud).getUpdatedMean(data);
            //         double test_var = getCloud(chosenCloud).getUpdatedVarianceK3(data, test_mean);
            //         // update the cloud stats
            //         getCloud(chosenCloud).updateStats(test_mean, test_var);
            //         // add the point to the cloud
            //         getCloud(chosenCloud).addPoint(index_point);
                
            //     }
            // }
            
            if (nothingCloud && isOutlier){
                // get the chosen cloud from the alfa list with the max value
                int maxIndex = std::distance(alfa.begin(), std::max_element(alfa.begin(), alfa.end()));

                // get chosen cloud
                DataCloud& c = getCloud(maxIndex);

                c.mean = data;

                c.adjust_variance_with_force(data);

                // Serial.println("Clouds reinitialized");

                c.addPoint(index_point);
            } else if (chosenCloud == -1 && isOutlier){
                
                //get the cloud with chosenCloud
                DataCloud& c = getCloud(chosenCloud);

                // update the cloud stats
                c.mean = data;

                c.adjust_variance_with_force(data);

                // add the point to the cloud
                c.addPoint(index_point);
            }

            return chosenCloud;
        }
    }

public:
    int index_point;
    // std::vector<DataCloud> clouds(3);
    // DataCloud clouds[3];  // Usando um array fixo
    
    // define a vector to store Dataclouds but wihtout initialize them
    std::vector<DataCloud> clouds;
    std::vector<std::vector<double>> classIndex;
    std::vector<int> argMax;
    std::vector<double> alfa = {0.0, 0.0};
    std::vector<double> relevanceList = {0.0, 0.0};
};