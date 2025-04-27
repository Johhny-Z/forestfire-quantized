#ifndef FORESTFIRE_HPP_
#define FORESTFIRE_HPP_

#include <cmath>
#include <limits>
#include <nlohmann/json.hpp>

#include <cadmium/modeling/celldevs/grid/cell.hpp>
#include <cadmium/modeling/celldevs/grid/config.hpp>
#include "forestfireState.hpp"

using namespace cadmium::celldevs;

class forestfire : public GridCell<forestfireState, double> {
    double wind_dir;
    double wind_spd;
    int myRow;
    int myCol;

public:
    forestfire(const std::vector<int>& id, 
            const std::shared_ptr<const GridCellConfig<forestfireState, double>>& config
          ): GridCell<forestfireState, double>(id, config) {
            config->rawCellConfig.at("wind").get_to(wind_dir);
            config->rawCellConfig.at("windspeed").get_to(wind_spd);
            myRow = id[0];
            myCol = id[1];
          }
    
    [[nodiscard]] forestfireState localComputation(
        forestfireState state,
        const std::unordered_map<std::vector<int>, NeighborData<forestfireState, double>>& neighborhood
    ) const override {
       
        double old_temp = state.temp;
        double new_temp = 0;
        double radians = (90 - wind_dir) * (M_PI / 180.0);
        double x = cos(radians);
        double y = sin(radians);

        for (const auto& [neighborId, neighborData] : neighborhood) {
            auto nState = neighborData.state;
            double v = neighborData.vicinity;
            double nTemp = nState->temp;
            
            int x_rel = neighborId[1] - myCol;
            int y_rel = neighborId[0] - myRow;
            
            double t_change = (1 + (wind_spd/30)*(x*x_rel - y*y_rel)) * nTemp * v;
            new_temp += t_change;
        }
        if (state.fire_status != 0) { 
            double delta = new_temp - old_temp;
            if (std::abs(delta) >= state.quantum) { 
                state.temp = new_temp; 
                state.sigma = 1.0;     
            } else {
                state.sigma = std::numeric_limits<double>::infinity(); 
            }
        }
        
        if (state.fire_status == 0) {
            state.temp = 300.0;
            state.sigma = std::numeric_limits<double>::infinity();
        }
        else if(state.fire_status == 1) {
            if(state.temp == 300.0){
				double unburned_temp = new_temp+0.213;
				state.temp = unburned_temp;
				if(state.temp == 300&&state.humidity>0.5){
					state.sigma = std::numeric_limits<double>::infinity();
				}
				else {
					double unburned_temp = new_temp+0.213;
                    state.last_temp = state.temp;
					state.temp = unburned_temp;
					state.fire_status = 2;	
					
				}
			}
            
        }    
        else if(state.fire_status == 2) {
            if(state.temp < 573.0) {
				double unburned_temp = new_temp+0.213;		
				state.temp = unburned_temp;
				if(state.temp == 300){
					state.sigma = std::numeric_limits<double>::infinity();
					state.fire_status = 1;	
				}	
				else {
					state.sigma = 1.0;
				}
			}	
			//if it ignites start the enthalpy formula, go to status 3
			else if(state.temp >= 573.0) {
				double burning_temp = new_temp+0.213+2.74*exp(-0.0019*(state.t_ig));
                state.last_temp = state.temp;		
				state.fire_status = 3;
				state.t_ig = 1.0;
				state.temp = burning_temp;
				
                
			}
        }
        else if(state.fire_status == 3) {
            if(state.temp > 333.0) {
                state.last_temp = state.temp;
                state.temp = new_temp + 0.213 + 2.74*exp(-0.0019*state.t_ig);
                state.t_ig += 1.0;
                
            }    
            else if((state.temp <= 333.0)||(state.humidity>0.6)) {
                state.last_temp = state.temp;
                state.temp = new_temp + 0.213;
                state.fire_status = 4;    
                
            }
        }
        else if(state.fire_status == 4) {
			double burned_temp = new_temp+0.213;		
			state.temp = burned_temp;
			if(state.temp == 300){
				state.sigma = std::numeric_limits<double>::infinity();
			}
			else {
				state.sigma = 1.0;
			}
		}

        // 更新量子化记录（新增）
        state.last_temp = state.temp;
        
        return state;
    }

    [[nodiscard]] double outputDelay(const forestfireState& state) const override {
        return (state.sigma == std::numeric_limits<double>::infinity()) ? 
            1e10 :  
            state.sigma;
    }
};

#endif