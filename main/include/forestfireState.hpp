#ifndef FORESTFIRE_STATE_HPP_
#define FORESTFIRE_STATE_HPP_

#include <iostream>
#include <nlohmann/json.hpp>


struct forestfireState {
	double temp; //temperature of the fire
	int fire_status; //0 - border, 1- inactive, 2 - not burned, 3 - burning, 4 - burned
	double t_ig; //ignition time
	double sigma;
	double ma;
	double last_temp;
	double quantum;
	double humidity;    
     	

	
	//! Default constructor function.
	forestfireState() : temp(300.0), fire_status(1), t_ig(0.0), sigma(1.0),ma(1.0),last_temp(300.0),quantum(0.03),humidity(0.2){}
};

//! It prints the temperature of the cell in an output stream.
std::ostream& operator<<(std::ostream& os, const forestfireState& x) {
	os << "<" << x.temp << ">";
	return os;
}
//! The simulator must be able to compare the equality of two state objects
bool operator!=(const forestfireState& x, const forestfireState& y) {
	return x.temp != y.temp ;
}
void from_json(const nlohmann::json& j, forestfireState& state) {
	j.at("temp").get_to(state.temp);
	j.at("fire_status").get_to(state.fire_status);
	j.at("t_ig").get_to(state.t_ig);
	j.at("quantum").get_to(state.quantum);
	j.at("humidity").get_to(state.humidity);
	state.last_temp = state.temp;
	
}

#endif 
