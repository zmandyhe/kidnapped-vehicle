/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */
#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"
#include "particle_filter.h"
#define EPS 0.00001

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
   //Initialize all particles to 
   //first position (based on estimates of x, y, theta and their uncertainties
   //from GPS) and all weights to 1. 
  if(is_initialized){ return;}
  //set the number of particles
  num_particles = 500;
  std::default_random_engine gen;
  //double std_x, std_y, std_theta;  // Standard deviations for x, y, and theta
  // TODO: Set standard deviations for x, y, and theta
  
  //Add random Gaussian noise to each particle.
  // This line creates a normal (Gaussian) distribution for x,y,theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  //generate particles with normal distribution
  for (int i=0; i<num_particles; i++){
    Particle particle;
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;
    
    particles.push_back(particle);
  }  
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate){
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  //This line creates a normal (Gaussian) distribution for x,y,theta
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  std::default_random_engine gen;
  //predict x,y,theta according to the formula
  for (int i=0; i<num_particles;i++){
    double theta = particles[i].theta;
    if(fabs(yaw_rate) < EPS){
      particles[i].x += velocity * delta_t * cos(theta);
      particles[i].y += velocity * delta_t * sin(theta);
    } else {
      particles[i].x += velocity/yaw_rate * (sin(theta + yaw_rate * delta_t) - sin(theta));
      particles[i].y += velocity/yaw_rate * (cos(theta) - cos(theta + yaw_rate * delta_t));
      particles[i].theta += yaw_rate * delta_t;
    }
    
   // adding noises
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta +=  dist_theta(gen);
  }  
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,vector<LandmarkObs>& observations) const{
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  // the predicted is a particle "P" estimated position to all lank marks
  unsigned int numObservations = observations.size();
  unsigned int numPredictions = predicted.size();
  
  for (unsigned int i = 0; i< numObservations; i++){
    double minDistance = std::numeric_limits<double>::max();
    int index_map = -1;
    for (unsigned int j=0; j<numPredictions; j++){
      double x_distance = observations[i].x - predicted[j].x;
      double y_distance = observations[i].y - predicted[j].y;
      double distance = x_distance * x_distance + y_distance * y_distance;
      if (distance < minDistance){
        minDistance = distance;
        index_map =  predicted[j].id;
      }
    }
    observations[i].id = index_map;
  }
}
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], const vector<LandmarkObs> &observations, 
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
  
  // transform each observation marker from the vehicle's coordinates to 
  //the map's coordinates, with respect to our particle using homogenous transformation matrix
  // transform to map x coordinate
  double stdLandmarkRange = std_landmark[0];
  double stdLandmarkBearing = std_landmark[1];
  for (int i = 0; i<num_particles; i++){
    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;
    //find landmarks from map data within paritcle's range
    double sensor_range_2 = sensor_range * sensor_range;
    vector<LandmarkObs> inRangeLandmarks;
    for (unsigned int j=0; j < map_landmarks.landmark_list.size();j++){
      float landmarkX = map_landmarks.landmark_list[j].x_f;
      float landmarkY = map_landmarks.landmark_list[j].y_f;
      int id = map_landmarks.landmark_list[j].id_i;
      double dX = x - landmarkX;
      double dY = y - landmarkY;
      if ( dX*dX + dY*dY <= sensor_range_2 ) {
        inRangeLandmarks.push_back(LandmarkObs{ id, landmarkX, landmarkY });
      }
//       if(dist(dX,dY,x,y) <= sensor_range){
//         inRangeLandmarks.push_back(LandmarkObs{id,landmarkX,landmarkY});
//       }
    }
    //transform observation coordinates
    vector<LandmarkObs> mappedObservations;
    for(unsigned int j=0; j<observations.size();j++){
      double xx = cos(theta)*observations[j].x - sin(theta)*observations[j].y + x;
      double yy = sin(theta)*observations[j].x + cos(theta)*observations[j].y + y;
      mappedObservations.push_back(LandmarkObs{observations[j].id, xx, yy});
    }
    //observation association to landmark in map
    ParticleFilter::dataAssociation(inRangeLandmarks,mappedObservations);
    //calculate weights
    particles[i].weight = 1.0;
    for (unsigned int j=0; j<mappedObservations.size();j++){
      double observationX = mappedObservations[j].x;
      double observationY = mappedObservations[j].y;
      
      int landmarkId = mappedObservations[j].id;
      
      double landmarkX, landmarkY;
      unsigned int k = 0;
      unsigned int nLandmarks = inRangeLandmarks.size();
      bool found = false;
      while(!found && k<nLandmarks){
        if(inRangeLandmarks[k].id == landmarkId){
          found = true;
          landmarkX = inRangeLandmarks[k].x;
          landmarkY = inRangeLandmarks[k].y;
        }
        k++;
      }
      //calculate weight
      double dX = observationX - landmarkX;
      double dY = observationY - landmarkY;
      double weight = ( 1/(2*M_PI*stdLandmarkRange*stdLandmarkBearing)) * exp((-dX*dX/(2*stdLandmarkRange*stdLandmarkRange)) + (-dY*dY/(2*stdLandmarkBearing*stdLandmarkBearing)));
//       double weight = ( 1/(2*M_PI*stdLandmarkRange*stdLandmarkBearing)) * exp((-dX*dX/(2*stdLandmarkRange*stdLandmarkRange) + (-dY*dY/(2*stdLandmarkBearing*stdLandmarkBearing))));
      //if(weight==0){
      if(weight < EPS){
        //particles[i].weight *= EPS;
        particles[i].weight = EPS;
      } else{
        particles[i].weight *= weight;
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
  vector<double> weights;
  double maxWeight = std::numeric_limits<double>::min();
  for(int i=0;i<num_particles;i++){
    weights.push_back(particles[i].weight);
    if(particles[i].weight > maxWeight){
      maxWeight = particles[i].weight;
    }
  }
  //creating distributions
  std::uniform_real_distribution<double> distDouble(0.0,maxWeight);
  std::uniform_int_distribution<int> distInt(0,num_particles-1);
  int index = distInt(gen);
  double b = distDouble(gen);
  double beta = 0.0;
  std::vector<Particle> resampledParticles;
  for(int i=0; i<num_particles; i++){
    beta += b * 2.0;
    while(beta > weights[index]){
      beta -= weights[index];
      index = (index+1) % num_particles;
    }
    resampledParticles.push_back(particles[index]);
  }
  particles = resampledParticles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                       const std::vector<double>& sense_x, 
                       const std::vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
  return particle;
}

std::string ParticleFilter::getAssociations(Particle best){
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  std::string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

std::string ParticleFilter::getSenseCoord(Particle best, std::string coord){
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