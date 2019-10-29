/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <cstddef>

#include "helper_functions.h"

#define EPS 0.00001
std::default_random_engine gen;
using namespace std;

//using std::string;
//using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
 
   if (is_initialized) {
    return;
   }
  
   n_particles = 169; // Initializing the number of particles
   // normal distribution of distribution x, y and theta (with std_x, std_x and std_theta separately)
   normal_distribution<double> dist_x(x, std[0]);
   normal_distribution<double> dist_y(y, std[1]);
   normal_distribution<double> angle_theta(theta, std[2]);
  
   // Generate particles with normal distribution with mean on GPS values.
   // for (int i = 0; i < n_particles; ++i) {
   for (int i = 0; i < n_particles; i++) {  
      //Using struct to make a particle structure and assign every information about each particles
      Particle particle;
      particle.id = i;
      particle.x = dist_x(gen);
      particle.y = dist_y(gen);
      particle.theta = angle_theta(gen);

      // assign weight=1 to each particle 
      particle.weight = 1.0;

      // add particle to ParticleFilter class =>  std::vector<Particle> particles;
      // with this method, every particle and vecotr particles can be generated.
      // add structure into a vector.
      particles.push_back(particle);
   }

  //after initialized, is_initialized should be true. If not, paricle fitler will always become initialized and uselessful.
  is_initialized = true;
 
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

   // normal distribution of distribution x with zero mean and each std.
   normal_distribution<double> dist_x(0, std_pos[0]);
   // normal distribution of distribution y with std_y
   normal_distribution<double> dist_y(0, std_pos[1]);
   // normal distribution of distribution theta with std_theta
   normal_distribution<double> angle_theta(0, std_pos[2]);
   
   // for loop
   for (int i = 0; i < n_particles; i++) {
         //double theta = particles[i].theta;
         if (fabs(yaw_rate) >= EPS) {
              particles[i].x = particles[i].x + (velocity / yaw_rate)*(sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
              particles[i].y = particles[i].y + (velocity / yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
              particles[i].theta = particles[i].theta + yaw_rate * delta_t;
         }
         else {// theta doesn't change
              particles[i].x = particles[i].x + velocity * delta_t *cos(particles[i].theta);
              particles[i].y = particles[i].y + velocity * delta_t *sin(particles[i].theta);
         }

         //add noise  to each particle in particles.
         particles[i].x = particles[i].x + dist_x(gen);
         particles[i].y = particles[i].y + dist_y(gen);
         particles[i].theta = particles[i].theta + angle_theta(gen);
   }
  
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

   int n_observation = observations.size();
   int n_predictions = predicted.size();

   for (int i = 0; i < n_observation; i++) {
         //for each observation
         //initializing the min distance as really big number
         double min_dis = numeric_limits<double>::max();
         
         //initializing the found map that is not in map , this is made for return the nearset measurement around GT.
         int id_in_map = -75;

         //complexity is o(ij);
         for (int j = 0; j < n_predictions; j++) {
               //distance calculation with helper function
               double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

               // if distance is smaller than the distance, then save the id , then iterate all the predicted value
               // finally find the most nearest precited to GT value. 
               if (distance < min_dis) {
                    min_dis = distance;
                    id_in_map = predicted[j].id;
               }
         }

         //assign the observed measurement to this particular landmark.
         //for vehicle, it means, this observation is belong to this landmark.
         observations[i].id = id_in_map;
   }  
  
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

   //std_landmark provides uncertainty of measurement of landmark in x,y direction . 
   double stdLandmarkRange = std_landmark[0];
   double stdLandmarkBearing = std_landmark[1];
   for (int i = 0; i < n_particles; i++) {
         double x = particles[i].x;
         double y = particles[i].y;
         double theta = particles[i].theta;

         //find landmarks in vehicle sensing range
         double sensor_range_2 = sensor_range * sensor_range;
         vector<LandmarkObs> inRangeLandmarks;
     
         for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
              float landmarkX = map_landmarks.landmark_list[j].x_f;
              float landmarkY = map_landmarks.landmark_list[j].y_f;
              int id = map_landmarks.landmark_list[j].id_i;
              double dX = x - landmarkX;
              double dY = y - landmarkY;

              //in this step, in range is constructed. After this step, we only calculate the landmarks in the range. 
              if (dX*dX + dY * dY <= sensor_range_2) {
                   inRangeLandmarks.push_back(LandmarkObs{ id, landmarkX, landmarkY });
              }
         }

         // Transfrom observation coodinates from vehicle coordinate to map (global) coordinate.
         vector<LandmarkObs> mappedObservations;
         
         //Rotation
         for (std::size_t j = 0; j< observations.size(); j++) {
               double xx = x + cos(theta)*observations[j].x - sin(theta) * observations[j].y;
               double yy = y + sin(theta)*observations[j].x + cos(theta) * observations[j].y;
               //using struct defined in helperfunction.h LandmarkObs, to make a after transition and rotation transformed observation data.
               mappedObservations.push_back(LandmarkObs{ observations[j].id, xx, yy });
         }

         //Observation association with landmark
         dataAssociation(inRangeLandmarks, mappedObservations);

         //reset the weight to 1.0
         particles[i].weight = 1.0;

         //calculate the weights
         for (std::size_t j = 0; j < mappedObservations.size(); j++) {
               double observationX = mappedObservations[j].x;
               double observationY = mappedObservations[j].y;
               int landmarkId = mappedObservations[j].id;
               double landmarkX, landmarkY;
               int k = 0;
               int nLandmarks = inRangeLandmarks.size();
               bool found = false;
               while (!found && k < nLandmarks) {
                    if (inRangeLandmarks[k].id == landmarkId) {
                         found = true;
                         landmarkX = inRangeLandmarks[k].x;
                         landmarkY = inRangeLandmarks[k].y;
                    }
                    k++;
               }
         
               //calculating weight
               double dX = observationX - landmarkX;
               double dY = observationY - landmarkY;

               //Since we assume the correlation between x direction and y direction is not exist, then rho in wiki is zero.
               //weight update
               double weight = (1 / (2 * M_PI*stdLandmarkRange*stdLandmarkBearing)) * exp(-(dX*dX / (2 * stdLandmarkRange*stdLandmarkRange) + (dY*dY / (2 * stdLandmarkBearing*stdLandmarkBearing))));

               //if weight equal to zero. then multiply to the EPS. 
               if (weight == 0) {
                    particles[i].weight = particles[i].weight*EPS;
               }
         
               //if weight doesn't equal to zero, then weight should be multiply by i times. because it is multivariate define.
               else {
                     particles[i].weight = particles[i].weight * weight;
               }
         }
   }
  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

   // Get weights and max weight.
   vector<double> weights;
   double max_weight = numeric_limits<double>::min();
   for(int i = 0; i < n_particles; i++) {
        weights.push_back(particles[i].weight);
        if ( particles[i].weight > max_weight ) {
              max_weight = particles[i].weight;
        }
   }
 
   // Creating distributions.
   uniform_real_distribution<float> dist_float(0.0, max_weight);
   uniform_int_distribution<int> dist_int(0, n_particles - 1);
   // Generating index.
   int index = dist_int(gen);
   double beta = 0.0;
  
   // the weight wheel
   vector<Particle> resampled_particles;
   for(int i = 0; i < n_particles; i++) {
        beta += dist_float(gen) * 2.0;
        while( beta > weights[index]) {
             beta -= weights[index];
             index = (index + 1) % n_particles;
        }
        resampled_particles.push_back(particles[index]);
   }
   particles = resampled_particles;
  
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  // and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}