#include<iostream>
#include<random>
#include<vector>
#include<fstream>
#include<cstdlib>
#include<Python.h>
using namespace std;

random_device rd; // Seed with a real random value, if available
mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

class Particle {
	private : 
		vector<double> x; //Current position at each time step
		vector<double> dx; //Position step at each time step 
		double t = 0; 	//Current time
	public : 
		void Initialization();
		void Brownian_Integral(int N, double dt); // Generating Brownian motion from current time step to next time step
		double Print_x_singletimestep(int i);
		vector<double> Print_x();
};


//Defining member function for class Particle
void Particle::Initialization() { x.push_back(0);}

void Particle::Brownian_Integral(int N, double dt) {
	//Generating Brownian Motion
	double mean = 0, std = sqrt(dt);
	
	//Generating normal distribution
	normal_distribution<> d(mean, std);
	for(int j = 1; j<N; j++) {
		//Sampling only 1 random number for Brownian motion : 
		double dx_current = d(gen);
		dx.push_back(dx_current);
		
		// Ensure the vector `x` is not empty
	    if (x.empty()) {
	        throw runtime_error("Vector x is empty. Ensure Initialization() is called before Brownian_Bridge().");
	    }
		
		double x_current = x.back() + dx_current;
		x.push_back(x_current);
	}

	//Generating Ito integral of Brownian bridge from 0 to t
	vector<double> Brownian = x;
	vector<double> integral = x;
	for(int i = 1; i<x.size(); i++) {
		integral[i] = integral[i-1] + 1/(1-(i-1)*dt)*dx[i-1];
		x[i] = (1 - i*dt)*integral[i];
	}
}

double Particle::Print_x_singletimestep(int i) { return x[i]; }

vector<double> Particle::Print_x() { return x; }


//Calculating statistical quantity from particles data
vector<double> Calculate_Mean(vector<Particle> particles, int N) {
	vector<double> mean;	//Mean from all time
	
	for(int i = 0; i<N; i++) {
		double sum = 0;
		for(int j = 0; j<particles.size(); j++) {
			sum += particles[j].Print_x_singletimestep(i);
		}
		mean.push_back(sum/particles.size());
	}
	
	return mean;
}

vector<vector<double>> Calculate_Fluctuation (vector<Particle> particles, int N, vector<double> mean, int num_particle) {
	vector<vector<double>> fluctuation;
	
	for(int i = 0; i<N; i++) {
		fluctuation.push_back(vector<double>());
		for(int j = 0; j<num_particle; j++) {
			fluctuation[i].push_back(particles[j].Print_x_singletimestep(i) - mean[i]);
		}
	}
	
	return fluctuation;
}

vector<double> Calculate_Variance(vector<vector<double>> fluctuation, int N, int num_particle) {
	vector<double> variance;
	
	for(int i = 0; i<N; i++) {
		double sum_variance = 0;
		
		for(int j = 0; j<num_particle; j++) {
			sum_variance += fluctuation[i][j]*fluctuation[i][j];
		}
		variance.push_back(sum_variance/num_particle);
	}
	
	return variance;
}

vector<vector<double>> Calculate_Covariance(vector<vector<double>> fluctuation, int N, int num_particles) {
    // Initialize an N x N covariance matrix with zeros
    vector<vector<double>> covariance_matrix(N, vector<double>(N, 0.0));

    // Calculate covariance for each pair of time steps (i, j)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {  // Only compute for j <= i, use symmetry for j > i
            double covariance = 0.0;

            // Summing the product of fluctuations for all particles
            for (int k = 0; k < num_particles; k++) {
                covariance += fluctuation[i][k] * fluctuation[j][k];
            }

            covariance /= num_particles;  // Normalize by number of particles
            covariance_matrix[i][j] = covariance;

            // Use symmetry: Cov(t_i, t_j) = Cov(t_j, t_i)
            if (i != j) {
                covariance_matrix[j][i] = covariance;
            }
        }
    }

    return covariance_matrix;
}

vector<vector<double>> Calculate_CorrelationCoeff(vector<vector<double>> covariance, vector<double> variance, int N, int num_particles) {
	vector<vector<double>> correlation_coeff(N, vector<double>(N, 0.0));
	
	for(int i = 0 ; i<N; i++) {
		for(int j = 0; j<=i; j++) {
			double denom = sqrt(variance[i]*variance[j]);
			if(denom != 0) { correlation_coeff[i][j] = covariance[i][j]/denom; }			
			if(i != j) { correlation_coeff[j][i] = correlation_coeff[i][j]; }
		}
	}
	return correlation_coeff;
}



//Output CSV
void OutputCSV_1D(string name, double dt, vector<Particle> particles, int N) {
	ofstream ofs;
	ofs.open(name);
	
	string upper_CSV = "t";
	for(int i = 0; i<particles.size(); i++) {
		upper_CSV += ",x" + to_string(i+1);
	}
	upper_CSV += "\n";
	ofs << upper_CSV;
	
	for (int i = 0; i<N; i++) {
		ofs << i*dt ;
		for(int j = 0; j<particles.size(); j++) {
			ofs << "," << particles[j].Print_x_singletimestep(i) ;
		}
		ofs << "\n";
	}
	ofs.close();
}

void OutputCSV_StatQuan1D(string name, double dt, vector<double> stat_quantity, int N, string quantity) {
	ofstream ofs;
	ofs.open(name);
	
	string upper_CSV = "t," + quantity + "\n";
	ofs<<upper_CSV;
	
	for (int i = 0; i<N; i++) {
		ofs << i*dt <<","<<stat_quantity[i]<<"\n";
	}
	ofs.close();
}

void OutputCSV_StatQuan2D(string name, double dt, vector<vector<double>> stat_quantity, int N, string quantity) {
	ofstream ofs;
	ofs.open(name);
	
	ofs << "t,";
    for (int j = 0; j < N; j++) {
        ofs << quantity << std::to_string(j);
        if (j != N - 1) {  // Add comma except for the last column
            ofs << ",";
        }
    }
    ofs << "\n";
	
	for(int i = 0; i<N; i++) {
		ofs<<i*dt;
		for(int j = 0; j<N; j++) {
			ofs<<","<<stat_quantity[i][j];
		}
		ofs<<"\n";
	}
	ofs.close();
}

//Plotting using Python
void call_Python() {
	// Python code for plotting
    const char* python_code = R"(
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.widgets import Slider

# Read the CSV file
data_particle = pd.read_csv('1D trajectory of particles (Brownian Bridge).csv')
data_mean = pd.read_csv('1D trajectory of particles (Brownian Bridge) Mean.csv')
data_variance = pd.read_csv('1D trajectory of particles (Brownian Bridge) Variance.csv')
data_covariance = pd.read_csv('1D trajectory of particles (Brownian Bridge) Covariance.csv')
data_correlation = pd.read_csv('1D trajectory of particles (Brownian Bridge) Correlation.csv')

# Time column
time_particle = data_particle['t']
time_mean = data_mean['t']
time_variance = data_variance['t']
time_covariance = data_covariance['t']
time_correlation = data_correlation['t']

#Especially data for covariance and correlation 
covariance = data_covariance.drop(columns=['t'])
correlation = data_correlation.drop(columns=['t'])

# Plot trajectories for each particle
plt.figure(1)
for column in data_particle.columns[1:]:
    plt.plot(time_particle, data_particle[column], label=column)

# Add labels, legend, and title
plt.xlabel('Time')
plt.ylabel('Position')
plt.title('Brownian Bridge 1D Trajectories')
plt.grid()
plt.show(block=False)  

x = np.zeros(data_mean.shape[0])
plt.figure(2)
plt.plot(time_mean, data_mean['mean'])
plt.plot(time_mean, x, 'k--', linewidth=2,)
plt.xlabel('Time')
plt.ylabel('Mean')
plt.title('Brownian Bridge 1D Trajectories (Mean)')
plt.grid()
plt.show(block=False)  

plt.figure(3)
plt.plot(time_variance, data_variance['variance'])
plt.plot(time_variance,  time_variance * (1-time_variance), 'k--', linewidth=2,)
plt.xlabel('Time')
plt.ylabel('Variance')
plt.title('Brownian Bridge 1D Trajectories (Variance)')
plt.grid()
plt.show(block=False)  

#For covariance
plt.figure(figsize=(8, 6))
plt.imshow(covariance, cmap='hot', interpolation='nearest')
plt.colorbar(label='Value')
plt.title('Heatmap of Covariance Matrix for Brownian Bridge')
plt.xlabel("Time")
plt.ylabel("Time")
plt.tight_layout()
plt.show(block=False)

#For correlation 
plt.figure(figsize=(8, 6))
plt.imshow(correlation, cmap='hot', interpolation='nearest')
plt.colorbar(label='Value')
plt.title('Heatmap of Correlation Matrix for Brownian Bridge')
plt.xlabel("Time")
plt.ylabel("Time")
plt.tight_layout()
plt.show(block=False)

#More on covariance and correlation matrix 
# Define the analytical covariance for Brownian Bridge
def analytical_covariance(s, t, T=1.0):
    # Convert inputs to NumPy arrays
    s = np.array(s)
    t = np.array(t)

    # Compute the covariance matrix
    min_st = np.minimum.outer(s, t)
    max_st = np.maximum.outer(s, t)
    covariance = min_st*(1-max_st)
    return covariance

# Analytical correlation matrix for Brownian Bridge
def analytical_correlation(s, t, T=1.0):
    """
    Computes the analytical correlation matrix for a Brownian Integral.
    """
    analytical_cov = analytical_covariance(s, t, T)
    stddev = np.sqrt(np.diag(analytical_cov))  # Standard deviations
    corr_matrix = analytical_cov / np.outer(stddev, stddev)
    corr_matrix[np.isnan(corr_matrix)] = 0  # Handle division by zero
    return corr_matrix

# Load numerical data from CSV
data_covariance = pd.read_csv('1D trajectory of particles (Brownian Bridge) Covariance.csv')
data_correlation = pd.read_csv('1D trajectory of particles (Brownian Bridge) Correlation.csv')

# Extract time points and numerical covariance/correlation matrices
time_points = data_covariance['t'].values
covariance = data_covariance.drop(columns=['t']).values
correlation = data_correlation.drop(columns=['t']).values

# Compute analytical covariance and correlation matrices
analytical_cov = analytical_covariance(time_points, time_points)
analytical_corr = analytical_correlation(time_points, time_points)

# Create the figure and subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
plt.subplots_adjust(bottom=0.25)
slider_ax = plt.axes([0.2, 0.1, 0.65, 0.03], facecolor='lightgoldenrodyellow')

# Initial time index
time_index = 0

# Initial numerical and analytical data
numerical_cov = covariance[time_index, :]
numerical_corr = correlation[time_index, :]
analytical_cov_row = analytical_cov[time_index, :]
analytical_corr_row = analytical_corr[time_index, :]

# Plot initial Covariance
num_cov_plot, = ax1.plot(time_points, numerical_cov, label="Numerical Covariance")
ana_cov_plot, = ax1.plot(time_points, analytical_cov_row, 'k--', linewidth=2, label="Analytical Covariance")
ax1.set_xlabel('Time')
ax1.set_ylabel('Covariance')
ax1.set_title(f'Covariance Comparison at Time t={time_points[time_index]:.2f}')
ax1.legend()
ax1.grid()

# Plot initial Correlation
num_corr_plot, = ax2.plot(time_points, numerical_corr, label="Numerical Correlation")
ana_corr_plot, = ax2.plot(time_points, analytical_corr_row, 'k--', linewidth=2, label="Analytical Correlation")
ax2.set_xlabel('Time')
ax2.set_ylabel('Correlation Coefficient')
ax2.set_title(f'Correlation Comparison at Time t={time_points[time_index]:.2f}')
ax2.legend()
ax2.grid()

# Slider for updating plots
time_slider = Slider(slider_ax, 'Time Index', 0, len(time_points) - 1, valinit=time_index, valstep=1)

# Update function for the slider
def update(val):
    time_index = int(time_slider.val)
    
    # Update numerical and analytical data
    numerical_cov = covariance[time_index, :]
    numerical_corr = correlation[time_index, :]
    analytical_cov_row = analytical_cov[time_index, :]
    analytical_corr_row = analytical_corr[time_index, :]

    # Update Covariance plot
    num_cov_plot.set_ydata(numerical_cov)
    ana_cov_plot.set_ydata(analytical_cov_row)
    ax1.set_title(f'Covariance Comparison at Time t={time_points[time_index]:.2f}')
    ax1.relim()
    ax1.autoscale_view()

    # Update Correlation plot
    num_corr_plot.set_ydata(numerical_corr)
    ana_corr_plot.set_ydata(analytical_corr_row)
    ax2.set_title(f'Correlation Comparison at Time t={time_points[time_index]:.2f}')
    ax2.relim()
    ax2.autoscale_view()

    fig.canvas.draw_idle()

# Attach the update function to the slider
time_slider.on_changed(update)

# Display the plots
plt.show(block=False)


input("Press Enter to close all plots...")
    )";

    // Execute the Python code
    PyRun_SimpleString(python_code);
}

int main() {	
	system("CLS");
	//Parameters
	double t_initial = 0, t_final = 1; 
	int N;
	int t_out = 1000;
	
	cout<<"From initial time = "<<t_initial<<" to final time = "<<t_final<<", how many partitions do you want? "; cin>>N;
	
	double dt = (t_final - t_initial)/(N-1);
	
	int num_particle_prompt;
	cout<<"Number of particles = "; cin>>num_particle_prompt;
    const int num_particle = num_particle_prompt;
    vector<Particle> particles(num_particle);
    
    //Numerical simulation
    for(int i = 0; i<num_particle; i++) {
    	if(num_particle > 100) {
    		if(i%t_out == 0) {
    			cout<<"Particle "<<i + 1<<endl;
    	
		    	//Initialization
		    	cout<<"Initialization : ";
		    	particles[i].Initialization();
		    	cout<<"Done\n";
		    	
		    	//Brownian motion
		    	cout<<"Simulating : ";
		    	particles[i].Brownian_Integral(N, dt);
				cout<<"Done\n\n";
			}
			else {  	
		    	//Initialization
		    	particles[i].Initialization();
		    	
		    	//Brownian motion
		    	particles[i].Brownian_Integral(N, dt);	
			}
    		
		}
		else {	
			cout<<"Particle "<<i + 1<<endl;
    	
	    	//Initialization
	    	cout<<"Initialization : ";
	    	particles[i].Initialization();
	    	cout<<"Done\n";
	    	
	    	//Brownian motion
	    	cout<<"Simulating : ";
	    	particles[i].Brownian_Integral(N, dt);
			cout<<"Done\n\n";
		}
	}
	
	//Calculate statistical quantity
	cout<<"Calculating Statistical Quantities : \n";
	
	cout<<"- Mean : ";
	vector<double> mean = Calculate_Mean(particles, N); cout<<"Done\n";
	
	cout<<"- Fluctuation : ";
	vector<vector<double>> fluctuation = Calculate_Fluctuation (particles, N, mean, particles.size()); cout<<"Done\n";
	
	cout<<"- Variance : ";
	vector<double> variance = Calculate_Variance(fluctuation, N, particles.size()); cout<<"Done\n";
	
	cout<<"- Covariance : ";
	vector<vector<double>> covariance = Calculate_Covariance(fluctuation, N, particles.size()); cout<<"Done\n";
	
	cout<<"- Correlation : ";
	vector<vector<double>> correlation = Calculate_CorrelationCoeff(covariance, variance, N, particles.size());
	cout<<"Done\n\n";
	
    //Output CSV
	cout<<"Outputing to CSV : \n";
	cout<<"- Trajectories : ";
	string name = "1D trajectory of particles (Brownian Bridge).csv"; OutputCSV_1D(name, dt, particles, N);   
	cout<<"Done\n"; 	
	
	cout<<"- Mean : ";
	name = "1D trajectory of particles (Brownian Bridge) Mean.csv"; OutputCSV_StatQuan1D(name, dt, mean, N, "mean");
	cout<<"Done\n";
	
	cout<<"- Variance : ";
	name = "1D trajectory of particles (Brownian Bridge) Variance.csv"; OutputCSV_StatQuan1D(name, dt, variance, N, "variance");
	cout<<"Done\n";
	
	cout<<"- Covariance : ";
	name = "1D trajectory of particles (Brownian Bridge) Covariance.csv"; OutputCSV_StatQuan2D(name, dt, covariance, N, "Cov");
	cout<<"Done\n";
	
	cout<<"- Correlation : ";
	name = "1D trajectory of particles (Brownian Bridge) Correlation.csv"; OutputCSV_StatQuan2D(name, dt, correlation, N, "Corr");
	
	cout<<"Done\n\n";
	
	//Calling Python for pkotting
	cout<<"Plotting using Python\n";
    Py_Initialize();
    call_Python();
    Py_Finalize();
    
    /*
    Type this in CMD
    g++ "Brownian Bridge.cpp" -IC:\Python313\include -LC:\Python313\libs -lpython313 -o "Brownian Bridge" -Wl,--enable-auto-import
	"Brownian Bridge.exe"
    */
}