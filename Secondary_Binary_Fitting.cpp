#include <iostream>
#include <istream>
#include <fstream>
#include <list>
#include <vector>
#include <cmath>

#include "xtensor/xtensor.hpp"
#include "xtensor/xarray.hpp"
#include "xtensor/xcsv.hpp"
#include "xtensor/xnpy.hpp"
#include "xtensor/xview.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xfixed.hpp"
#include "xtensor/xmath.hpp" 
#include "xsimd/xsimd.hpp"
#include "xtensor/xstrides.hpp"
#include "xtensor/xsort.hpp"
#include "xtensor/xadapt.hpp"
using namespace std;

vector<vector<double>> BinaryBase(vector<int> primary_coordinates, vector<int> secondary_coordinates) 
{
	vector<int> w1(5);
	vector<int> w2(5);
	
	w1(primary_coordinates(0) - 1, primary_coordinates(1) - 1) = 1;
	w1(primary_coordinates(0) - 1, primary_coordinates(1)) = 1;
	w1(primary_coordinates(0) - 1, primary_coordinates(1) + 1) = 1;
	w1(primary_coordinates(0), primary_coordinates(1) - 1) = 1;
	w1(primary_coordinates(0), primary_coordinates(1)) = 1;
	w1(primary_coordinates(0), primary_coordinates(1) + 1) = 1;
	w1(primary_coordinates(0) + 1, primary_coordinates(1) - 1) = 1;
	w1(primary_coordinates(0) + 1, primary_coordinates(1)) = 1;
	w1(primary_coordinates(0) + 1, primary_coordinates(1) + 1) = 1;

	w2(secondary_coordinates(0) - 1, secondary_coordinates(1) - 1) = 1;
	w2(secondary_coordinates(0) - 1, secondary_coordinates(1)) = 1;
	w2(secondary_coordinates(0) - 1, secondary_coordinates(1) + 1) = 1;
	w2(secondary_coordinates(0), secondary_coordinates(1) - 1) = 1;
	w2(secondary_coordinates(0), secondary_coordinates(1)) = 1;
	w2(secondary_coordinates(0), secondary_coordinates(1) + 1) = 1;
	w2(secondary_coordinates(0) + 1, secondary_coordinates(1) - 1) = 1;
	w2(secondary_coordinates(0) + 1, secondary_coordinates(1)) = 1;
	w2(secondary_coordinates(0) + 1, secondary_coordinates(1) + 1) = 1;
	
	auto w = w1 || w2;
	vector<vector<int>> w_array = {w, w1, w2};
	
	return w_array;
}

vector<int> BinaryPSF(vector<vector<double>> array_of_PSFs, vector<int> secondary_coordinates, vector<double> image_psf, double flux, vector<int> w)
{
    cout << secondary_coordinates << endl;
	vector<double> scale_array;
    for(i= 0.95; i > 0; i = i-0.05;){
        scale_array.push_back(i);
    }
	int models_length = array_of_PSFs.size()[0];
	int scale_length = scale_array.size()[0];
	int chi_array_length = models_length * models_length * scale_length;
	
	vector<double> chi_array;
	chi_array.reserve(chi_array_length);
	
	int cx = (array_of_PSFs(0).size()[1] - 1) * .5;
	int cy = (array_of_PSFs(0).size()[0] - 1) * .5;
	cout << "line 66" << endl;
	for (int i = 0; i < models_length; ++i) { 
		vector<double> model_PSF = array_of_PSFs(i);
		vector<double> PSF1 = (xt::view(model_PSF, xt::range(cy - 2, cy + 3 ), xt::range(cx - 2, cx + 3)));
		for (int j = 0; j < models_length; ++j) {
			vector<double> model_PSF = array_of_PSFs(j);
			int rangeY_1 = cy - secondary_coordinates(0);
			int rangeY_2 = cy - secondary_coordinates(0) + 5;
			int rangeX_1 = cx - secondary_coordinates(1);
			int rangeX_2 = cx - secondary_coordinates(1) + 5;
			
			vector<double> PSF2 = (xt::view(model_PSF, xt::range(rangeY_1, rangeY_2), xt::range(rangeX_1, rangeX_2)));	
			for (int k = 0; k < scale_length; ++k) {
				double relative_flux = scale_array(k);
				vector<double> a = PSF1 * relative_flux;
				vector<double> b = PSF2 * (1 - relative_flux);
				vector<double> temp_array = (a + b) * flux;
				vector<double> chi_PSF = w * (image_psf - temp_array);
				double chi_value = sum(chi_PSF * chi_PSF)(0);
				chi_array.push_back(chi_value);	
			}
		}	
	} 
    cout << "line 89" << endl;
	auto chi = xt::adapt(chi_array, {100, 100, 19});
	vector<int> best_PSFs = xt::adapt(xt::unravel_index(xt::argmin(chi)(0), chi.size(), xt::layout_type::row_major), {3});
	cout << "line 91: " << best_PSFs << endl;
	return best_PSFs;
}

vector<double> BinaryRelativeFlux(vector<vector<double>> array_of_PSFs, vector<int> secondary_coordinates, vector<double> image_psf, vector<double> flux, vector<int> w, int best_primary, int best_secondary) 
{
	//ref = relative flux
	vector<int> primary_ref(9);
	vector<int> secondary_ref(9);
	
	if (best_primary < 10) {
		best_primary += 10;
	}
    	if (best_primary >= 90) {
        	best_primary -= 10;
        }
    	if (best_primary % 10 == 0) {
        	best_primary += 1;
	}
   	if (best_primary % 10 == 9) {
        	best_primary-= 1;
	}
	
	if (best_secondary < 10) {
		best_secondary  += 10;
	}
    	if (best_secondary >= 90) {
        	best_secondary -= 10;
        }
    	if (best_secondary % 10 == 0) {
        	best_secondary  += 1;
	}
   	if (best_secondary % 10 == 9) {
        	best_secondary  -= 1;
	}
	
	primary_ref(0) = best_primary - 11;
	primary_ref(1) = best_primary - 10;
	primary_ref(2) = best_primary - 9;
	primary_ref(3) = best_primary - 1;
	primary_ref(4) = best_primary;
	primary_ref(5) = best_primary + 1;
	primary_ref(6) = best_primary + 9;
	primary_ref(7) = best_primary + 10;
	primary_ref(8) = best_primary + 11;
	
	secondary_ref(0) = best_secondary - 11;
	secondary_ref(1) = best_secondary - 10;
	secondary_ref(2) = best_secondary - 9;
	secondary_ref(3) = best_secondary - 1;
	secondary_ref(4) = best_secondary;
	secondary_ref(5) = best_secondary + 1;
	secondary_ref(6) = best_secondary + 9;
	secondary_ref(7) = best_secondary + 10;
	secondary_ref(8) = best_secondary + 11;
	
	vector<double> scale_array = xt::arange<double>(.99, 0, -.01);
	
	int primary_length = primary_ref.size()[0];
	int secondary_length = secondary_ref.size()[0];
	int scale_length = scale_array.size()[0];
	int chi_array_length = primary_length * secondary_length * scale_length;
	vector<double> chi_array;
	chi_array.reserve(chi_array_length);
	
	for (int i = 0; i < primary_length; ++i ) {
		int primary = primary_ref(i);
		vector<double> model_PSF = array_of_PSFs(primary);
		int cx = (model_PSF.size()[1] - 1) * .5;
		int cy = (model_PSF.size()[0] - 1) * .5;
		vector<double> PSF1 = (xt::view(model_PSF, xt::range(cy - 2, cy + 3 ), xt::range(cx - 2, cx + 3)));	
		for (int j = 0; j < secondary_length; ++j) {
			int secondary = secondary_ref(j);
			vector<double> model_PSF = array_of_PSFs(secondary);
			int rangeY_1 = cy - secondary_coordinates(0);
			int rangeY_2 = cy - secondary_coordinates(0) + 5;
			int rangeX_1 = cx - secondary_coordinates(1);
			int rangeX_2 = cx - secondary_coordinates(1) + 5;
			
			vector<double> PSF2 = (xt::view(model_PSF, xt::range(rangeY_1, rangeY_2), xt::range(rangeX_1, rangeX_2)));	
			for (int k = 0; k < scale_length; ++k) {
				double relative_flux = scale_array(k);
				vector<double> a = PSF1 * relative_flux;
				vector<double> b = PSF2 * (1 - relative_flux);
				vector<double> temp_array = (a + b) * flux;
				vector<double> chi_PSF = w * (image_psf - temp_array);
				double chi_value = xt::sum(chi_PSF * chi_PSF, xt::evaluation_strategy::immediate)(0);
				chi_array.push_back(chi_value);
			}
		}
	}
	
	vector<double> chi = xt::adapt(chi_array, {9, 9, 99});
	vector<int> best_PSFs = xt::adapt(xt::unravel_index(xt::argmin(chi)(0), chi.size(), xt::layout_type::row_major), {3});
	double best_primary_ref = primary_ref(best_PSFs(0));
	double best_secondary_ref = secondary_ref(best_PSFs(1));
	double best_ref = scale_array(best_PSFs(2));
	vector<double> final_best_PSF = {best_primary_ref, best_secondary_ref, best_ref};
	
	return final_best_PSF;
}

double BinaryFit(vector<double> model_primary, vector<double> model_secondary, vector<int> secondary_coordinates, vector<double> image_psf, double base, vector<int> w, double relative_flux) 
{
	int cx = (model_primary.size()[1] - 1) * .5;
	int cy = (model_primary.size()[0] - 1) * .5;
	int rangeY_1 = cy - secondary_coordinates(0);
	int rangeY_2 = cy - secondary_coordinates(0) + 5;
	int rangeX_1 = cx - secondary_coordinates(1);
	int rangeX_2 = cx - secondary_coordinates(1) + 5;
	
	vector<double> PSF1 = (xt::view(model_primary, xt::range(cy - 2, cy + 3 ), xt::range(cx - 2, cx + 3)));		
	vector<double> PSF2 = (xt::view(model_secondary, xt::range(rangeY_1, rangeY_2), xt::range(rangeX_1, rangeX_2)));
	vector<double> flux_array;
	vector<double> chi_array;
	
	//change "range" to something like "Scale"
	vector<double> scale_array = xt::arange<double>(0, 10.01, 0.01);
	int scale_length = scale_array.size()[0];
	flux_array.reserve(scale_length);
	chi_array.reserve(scale_length);
	vector<double> PSF_sum = PSF1 * relative_flux + PSF2 * (1 - relative_flux);
	double scale;
	for (int i = 0; i < scale_length; ++i) {
		scale = scale_array(i);
		vector<double> scale_PSF = PSF_sum * scale * base;
		flux_array.push_back(scale * base);
		vector<double> chi_PSF = w * (image_psf - scale_PSF);
		double chi_value = xt::sum(chi_PSF * chi_PSF, xt::evaluation_strategy::immediate)(0);
		chi_array.push_back(chi_value);
	}
		
	auto chi_squareds = xt::adapt(chi_array, {scale_length});
	auto fluxes = xt::adapt(flux_array, {scale_length});
	double min_chi = xt::amin(chi_squareds)(0);
	int min_chi_index = 0;
	int chi_length = chi_squareds.size()[0];
	//USE A VARIABLE FOR CHI_SQUAREDS.SHAPE, CHECK IN THE OTHER FILES TOO
	for (int i = 0; i < chi_length; ++i) {
		if (chi_squareds(i) == min_chi) {
			min_chi_index = i;
			break;
		}
	}	
	
	double best_flux = fluxes(min_chi_index);
	
	return best_flux;	 
}

double SecondaryBinaryFit(vector<double> primary, vector<double> secondary, vector<int> secondary_coordinates, vector<int >w2, vector<double> image_psf, double flux, double relative_flux)
{
	int cx = (primary.size()[1] - 1) * .5;
	int cy = (primary.size()[0] - 1) * .5;
	int rangeY_1 = cy - secondary_coordinates(0);
	int rangeY_2 = cy - secondary_coordinates(0) + 5;
	int rangeX_1 = cx - secondary_coordinates(1);
	int rangeX_2 = cx - secondary_coordinates(1) + 5;
	
	vector<double> PSF1 = (xt::view(primary, xt::range(cy - 2, cy + 3), xt::range(cx - 2, cx + 3)));
	vector<double> PSF2 = (xt::view(secondary, xt::range(rangeY_1, rangeY_2), xt::range(rangeX_1, rangeX_2)));
	vector<double> scale_array = xt::arange<double>(0, 10.0001, .0001);
	int scale_length = scale_array.size()[0];
	vector<double> flux_array;
	vector<double> chi_array;
	flux_array.reserve(scale_length);
	chi_array.reserve(scale_length);
	
	vector<double> relative_PSF2 = w2 * (image_psf - PSF1 * flux * relative_flux);
	for (int i = 0; i < scale_length; ++i)
	{
		double scale = scale_array(i);
		vector<double> scale_PSF2 =  w2 * PSF2 * scale * flux;
		flux_array.push_back(scale * flux);
		vector<double> chi_PSF = w2 * (relative_PSF2 - scale_PSF2);
		double chi_value = xt::sum(chi_PSF * chi_PSF, xt::evaluation_strategy::immediate)(0);
		chi_array.push_back(chi_value);
	}

	auto chi_squareds = xt::adapt(chi_array, {scale_length});
	auto fluxes = xt::adapt(flux_array, {scale_length});
	double min_chi = xt::amin(chi_squareds)(0);
	int min_chi_index = 0;
	int chi_length = chi_squareds.size()[0];
	for (int i = 0; i < chi_length ; ++i) {
		if (chi_squareds(i) == min_chi) {
			min_chi_index = i;
			break;
		}
	}	
	
	double best_flux = fluxes(min_chi_index);
	
	return best_flux;
}

double RelativeFluxes(vector<double> PSF_1, vector<int> coordinates_1, vector<double> PSF_2, vector<int> coordinates_2, vector<int> w, vector<double> image_psf, double base, int flux) 
{
	int cx = (PSF_1.size()[1] - 1) * .5;
	int cy = (PSF_2.size()[0] - 1) * .5;
	
	int P_rangeY_1 = cy - coordinates_1(0);
	int P_rangeY_2 = cy - coordinates_1(0) + 5;
	int P_rangeX_1 = cx - coordinates_1(1);
	int P_rangeX_2 = cx - coordinates_1(1) + 5;
	
	int S_rangeY_1 = cy - coordinates_2(0);
	int S_rangeY_2 = cy - coordinates_2(0) + 5;
	int S_rangeX_1 = cx - coordinates_2(1);
	int S_rangeX_2 = cx - coordinates_2(1) + 5;

	vector<double> PSF1 = (xt::view(PSF_1, xt::range(P_rangeY_1, P_rangeY_2), xt::range(P_rangeX_1, P_rangeX_2)));
	vector<double> PSF2 = (xt::view(PSF_2, xt::range(S_rangeY_1, S_rangeY_2), xt::range(S_rangeX_1, S_rangeX_2)));
	
	vector<double> scale_array = xt::arange<double>(0, 10.0001, .0001);
	int scale_length = scale_array.size()[0];
	vector<double> flux_array;
	vector<double> chi_array;
	flux_array.reserve(scale_length);
	chi_array.reserve(scale_length);
	
	auto rPSF = w * (image_psf - PSF2 * flux);
	for (int i = 0; i < scale_length; ++i) {
		double scale = scale_array(i);
		vector<double> scale_PSF = w * PSF1 * scale * base;
		flux_array.push_back(scale * base);
		vector<double> chi_PSF = w * (rPSF - scale_PSF);
		double chi_value = xt::sum(chi_PSF * chi_PSF, xt::evaluation_strategy::immediate)(0);
		chi_array.push_back(chi_value);
	}
	
	auto chi_squareds = xt::adapt(chi_array, {scale_length});
	auto fluxes = xt::adapt(flux_array, {scale_length});
	double min_chi = xt::amin(chi_squareds)(0);
	int min_chi_index = 0;
	int chi_length = chi_squareds.size()[0];
	for (int i = 0; i < chi_length ; ++i) {
		if (chi_squareds(i) == min_chi) {
			min_chi_index = i;
			break;
		}
	}	
	
	double best_flux = fluxes(min_chi_index);
	
	return best_flux;
}

vector<double> Read_File(string filename){
    vector<double> vector; 
    short loop=0; 
    double line; 
    ifstream myfile(filename); 
    if (myfile.is_open()) 
    {
        while (! myfile.eof() ) 
        {
            getline(myfile,line); 
            vector.push_back(line);
            cout << array[loop] << endl; 
            loop++;
        }
        myfile.close(); 
    }
    return vector;
}

vector<int> slicing(vector<int>& arr,
                    int X, int Y)
{
 
    // Starting and Ending iterators
    auto start = arr.begin() + X;
    auto end = arr.begin() + Y + 1;
 
    // To store the sliced vector
    vector<int> result(Y - X + 1);
 
    // Copy vector using copy function()
    copy(start, end, result.begin());
 
    // Return the final sliced vector
    return result;
}
//MAIN 
int main() 
{
	//read in the image and model data, and the parameters we set earlier
	vector<double> image_psf = Read_File("tempdata/image_psf.txt");   

	int array_length = 100;
	vector<vector<double>> array_of_PSFs(100);
//	vector<vector<double>> array_of_PSFs = xt::load_npy<double>("tempdata/array_of_PSFs.txt");
    	for (int i = 1; i <= array_length; ++i) {
    		string filename = "tempdata/PSFmodel_" + to_string(i) + ".txt";
    		array_of_PSFs(i - 1) = Read_File(filename);
    	}
    	
	auto coordinates = Read_File("tempdata/sec_coords.txt");
	auto center_load = Read_File("tempdata/prim_coords.txt");
	vector<int> center = {(int) center_load(0,0), (int) center_load(1,1)};
    	
    	vector<int> primary_coordinates = {2,2};
    	vector<int> secondary_coordinates = vector<double>(2);
    	
    	//arrays to store the output values
    	vector<double> residual_error_array(coordinates.size()[0]);
	vector<int> best_primary_array(coordinates.size()[0]);
	vector<double> flux_primary_array(coordinates.size()[0]);
	vector<int> best_secondary_array(coordinates.size()[0]);
	vector<double> flux_secondary_array(coordinates.size()[0]);
	vector<double> flux_sum_array(coordinates.size()[0]);
    	
    	//do the binary fit
    	cout << "Beginning Second Binary Fit" << endl;
        cout << "coordinates" << coordinates << "center" << center << endl;
    
    //The coordinates provided are the 9 pixels surrounding the centroid
    //For the second binary fit, we want to provide the 9 subpixels that surround the best fit primary and the secondary, in this case, 141.73,150.45 and 143.23,150.57
    	for (int i = 0; i < coordinates.size()[0]; ++i) {
            cout << "loop " << i << endl;
    		vector<int> secondary = coordinates[i];
    		secondary_coordinates[0] = secondary[0] - center[0] + 2;
    		secondary_coordinates[1] = secondary[1] - center[1] + 2;
    		//Create a new Binary Base from the new coordinates given
    		vector<vector<double>> w_array = BinaryBase(primary_coordinates, secondary_coordinates);
		vector<double> w = w_array[0];
		vector<int> w1 = w_array[1];
		vector<int> w2 = w_array[2];
		auto sqrt_weights = sqrt(abs(image_psf));
		
    		int iteration = 0;
    		int best_primary = 0;
    		int best_secondary = 0;
    		double best_RF = 0.0;
    		double base = sum(w * image_psf)(0);
    		
    		double comp = 0;
    		double flux = base * 1;
    		cout << "line 395" << endl;
    		while (abs(comp - flux) > .000001 && iteration < 10) {
    			comp = flux;
                cout << "line 398" << endl;
    			vector<int> best_PSFs =  BinaryPSF(array_of_PSFs, secondary_coordinates, image_psf, flux, sqrt_weights);
                cout << "line 401" << endl;
    			best_primary = best_PSFs(0);
    			best_secondary = best_PSFs(1);
                cout << "line 402" << endl;
    			vector<double> best_PSFs_ref = BinaryRelativeFlux(array_of_PSFs, secondary_coordinates, image_psf, flux, sqrt_weights, best_primary, best_secondary);
    			best_primary = best_PSFs_ref(0);
    			best_secondary = best_PSFs_ref(1);
    			best_RF = best_PSFs_ref(2);
    			cout << "line 407" << endl;
			flux = BinaryFit(array_of_PSFs[best_primary], array_of_PSFs[best_secondary], secondary_coordinates, image_psf, base, w, best_RF);
                cout << "line 409" << endl;
    			++iteration ;
    		}
    		
    		vector<double> primary_PSF = array_of_PSFs[best_primary];
    		vector<double> secondary_PSF = array_of_PSFs[best_secondary];

    		double comp_primary = 1.0;
    		double comp_secondary = 1.0;
    		double flux_primary = 2.0;

    		double flux_secondary = SecondaryBinaryFit(primary_PSF, secondary_PSF, secondary_coordinates, w2, image_psf, flux, best_RF);

    		iteration = 0;
    		
    		//Final WHILE LOOP 
		while (abs(comp_primary - flux_primary) > .000001 && abs(comp_secondary - flux_secondary) > .000001) {
			comp_primary = flux_primary;
			comp_secondary = flux_secondary;
			if (iteration >= 10) {
				break;
			}
			//CREATE rfpsf FUNCTION
			flux_primary = RelativeFluxes(primary_PSF, primary_coordinates, secondary_PSF, secondary_coordinates, w1, image_psf, base, flux_secondary);
			flux_secondary = RelativeFluxes(secondary_PSF, secondary_coordinates, primary_PSF, primary_coordinates, w2, image_psf, base, flux_primary);
			++iteration;
            
		}
		
		int cx = (primary_PSF.size()[1] - 1) * .5;
		int cy = (primary_PSF.size()[0] - 1) * .5;
		
		int P_rangeY_1 = cy - primary_coordinates[0];
		int P_rangeY_2 = cy - primary_coordinates[0] + 5;
		int P_rangeX_1 = cx - primary_coordinates[1];
		int P_rangeX_2 = cx - primary_coordinates[1] + 5;
		
		int S_rangeY_1 = cy - secondary_coordinates[0];
		int S_rangeY_2 = cy - secondary_coordinates[0] + 5;
		int S_rangeX_1 = cx - secondary_coordinates[1];
		int S_rangeX_2 = cx - secondary_coordinates[1] + 5;
       
        vector<vector<double>> PSF1;
        PSF1.reserve(P_rangeY_2-P_rangeY_1);
        vector<vector<double>> PSF2;
        PSF2.reserve(S_rangeY_2-S_rangeY_1);
        for (size_t i = P_range_Y_1; i < P_range_Y_2; ++i) {
            PSF1.emplace_back(primary_PSF[i].begin() + P_rangeX_1, primary_PSF[i].begin() + P_rangeX_2);
        }
        for (size_t i = S_range_Y_1; i < S_range_Y_2; ++i) {
            PSF2.emplace_back(primary_PSF[i].begin() + S_rangeX_1, primary_PSF[i].begin() + S_rangeX_2);
        }
		//vector<double> PSF1 = slicing(primary_PSF, xt::range(P_rangeY_1, P_rangeY_2), xt::range(P_rangeX_1, P_rangeX_2)));
		//vector<double> PSF2 = (xt::view(secondary_PSF, xt::range(S_rangeY_1, S_rangeY_2), xt::range(S_rangeX_1, S_rangeX_2)));

		vector<double> PSF_sum = PSF1 * flux_primary + PSF2 * flux_secondary;
		double pre_residual = abs(image_psf - PSF_sum);
		double residual = sum(w * pre_residual)[0];
		double residual_error = residual / sum(w)[0];
		
		cout << "psf1: " << best_primary << " flux1: " << flux_primary << " psf2: " << best_secondary << " flux2: " << flux_secondary << " error: " << residual_error << endl;
		//append each variable to a vector of that variable
		residual_error_array[i] = residual_error;
		best_primary_array[i] = best_primary;
		flux_primary_array[i] = flux_primary;
		best_secondary_array[i] = best_secondary;
		flux_secondary_array[i] = flux_secondary;
		flux_sum_array[i] = flux_primary + flux_secondary;
    	}
	
	cout << endl;
    	vector<int> center_array = {{center[0], center[1]}, {center[0], center[1]}, {center[0], center[1]}, {center[0], center[1]}, {center[0], center[1]}, {center[0], center[1]}, {center[0], center[1]}, {center[0], center[1]}, {center[0], center[1]}};
    	vector<double> secondary_array = coordinates;
    	
    	//send the "outputs" as seperate files since they are different data types and we can read them all in to python for the calculation of angles and separation - it only takes a second :)
    ofstream outFile1("tempdata/residual_error_array.txt");
    for (const auto &e : residual_error_array) outFile1 << e << "\n";
    ofstream outFile2("tempdata/center_array.txt");
    for (const auto &e : center_array) outFile2 << e << "\n";
    ofstream outFile3("tempdata/secondary_array.txt");
    for (const auto &e : secondary_array) outFile3 << e << "\n";
    ofstream outFile4("tempdata/best_primary_array.txt");
    for (const auto &e : best_primary_array) outFile4 << e << "\n";
    ofstream outFile5("tempdata/flux_primary_array.txt");
    for (const auto &e : flux_primary_array) outFile5 << e << "\n";
    ofstream outFile6("tempdata/best_secondary_array.txt");
    for (const auto &e : best_secondary_array) outFile6 << e << "\n";
    ofstream outFile7("tempdata/flux_sum_array.txt");
    for (const auto &e : flux_sum_array) outFile7 << e << "\n";
    ofstream outFile8("tempdata/flux_secondary_array.txt");
    for (const auto &e : flux_secondary_array) outFile8 << e << "\n";
    
	return 0;
}
