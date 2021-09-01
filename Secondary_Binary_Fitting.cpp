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

//A copy of the BinaryBase function in Binary_Fitting
//Returns a 9x9 vector of binary values
vector<vector<int>> BinaryBase(pair<int,int> primary_coordinates, pair<int,int> secondary_coordinates) 
{
	vector<vector<int>> w1 {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
	vector<vector<int>> w2 {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
    vector<vector<int>> w{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}}; 
    
    w1[primary_coordinates.first - 1][ primary_coordinates.second - 1] = 1;
	w1[primary_coordinates.first - 1][ primary_coordinates.second] = 1;
	w1[primary_coordinates.first - 1][ primary_coordinates.second + 1] = 1;
	w1[primary_coordinates.first][ primary_coordinates.second - 1] = 1;
	w1[primary_coordinates.first][ primary_coordinates.second] = 1;
	w1[primary_coordinates.first][ primary_coordinates.second + 1] = 1;
	w1[primary_coordinates.first+ 1][ primary_coordinates.second - 1] = 1;
	w1[primary_coordinates.first + 1][ primary_coordinates.second] = 1;
	w1[primary_coordinates.first + 1][ primary_coordinates.second + 1] = 1;

    w[primary_coordinates.first - 1][ primary_coordinates.second - 1] = 1;
	w[primary_coordinates.first - 1][ primary_coordinates.second] = 1;
	w[primary_coordinates.first - 1][ primary_coordinates.second + 1] = 1;
	w[primary_coordinates.first][ primary_coordinates.second - 1] = 1;
	w[primary_coordinates.first][ primary_coordinates.second] = 1;
	w[primary_coordinates.first][ primary_coordinates.second + 1] = 1;
	w[primary_coordinates.first + 1][ primary_coordinates.second - 1] = 1;
	w[primary_coordinates.first + 1][ primary_coordinates.second] = 1;
	w[primary_coordinates.first + 1][ primary_coordinates.second + 1] = 1;
    
	w2[secondary_coordinates.first - 1][ secondary_coordinates.second - 1] = 1;
	w2[secondary_coordinates.first - 1][ secondary_coordinates.second] = 1;
	w2[secondary_coordinates.first - 1][ secondary_coordinates.second + 1] = 1;
	w2[secondary_coordinates.first][ secondary_coordinates.second - 1] = 1;
	w2[secondary_coordinates.first][ secondary_coordinates.second] = 1;
	w2[secondary_coordinates.first][ secondary_coordinates.second + 1] = 1;
	w2[secondary_coordinates.first + 1][ secondary_coordinates.second - 1] = 1;
	w2[secondary_coordinates.first + 1][ secondary_coordinates.second] = 1;
	w2[secondary_coordinates.first + 1][ secondary_coordinates.second + 1] = 1;
    
	w[secondary_coordinates.first - 1][ secondary_coordinates.second - 1] = 1;
	w[secondary_coordinates.first - 1][ secondary_coordinates.second] = 1;
	w[secondary_coordinates.first - 1][ secondary_coordinates.second + 1] = 1;
	w[secondary_coordinates.first][ secondary_coordinates.second - 1] = 1;
	w[secondary_coordinates.first][ secondary_coordinates.second] = 1;
	w[secondary_coordinates.first][ secondary_coordinates.second + 1] = 1;
	w[secondary_coordinates.first + 1][ secondary_coordinates.second - 1] = 1;
	w[secondary_coordinates.first + 1][ secondary_coordinates.second] = 1;
	w[secondary_coordinates.first + 1][ secondary_coordinates.second + 1] = 1;
	

	vector<vector<int>> w_array{};
    w_array.insert( w_array.end(), w.begin(), w.end() );
    w_array.insert( w_array.end(), w1.begin(), w1.end() );
    w_array.insert( w_array.end(), w2.begin(), w2.end() );
	
	return w_array;
}

//This function Reshapes a 1D Vector into a 3D vector of the given shape parameters in the double form
vector<vector<vector<double>>> Reshape3D(vector<double> chi_array, int shape1, int shape2, int shape3){
    vector<vector<vector<double>>> output(shape1,vector<vector<double>>(shape2,vector<double>(shape3)));
    for(int i=0;i<shape1;i++){
        for(int j=0;j<shape2;j++){
            for(int k=0;k<shape3;k++){
                int num = i*shape3*shape2+shape3*j+k;
                cout << chi_array.size() << " " << num;
                output[i][j][k] = chi_array[num];
            }
        }
    }
    return output;
}

//This function Reshapes a 1D Vector into a 3D vector of the given shape parameters in the int form
vector<vector<vector<int>>> Reshape3DInt(vector<double> chi_array, int shape1, int shape2, int shape3){
    vector<vector<vector<int>>> output(shape1,vector<vector<int>>(shape2,vector<int>(shape3)));
    for(int i=0;i<shape1;i++){
        for(int j=0;j<shape2;j++){
            for(int k=0;k<shape3;k++){
                int num = i*shape3*shape2+shape3*j+k;
                output[i][j][k] = chi_array[num];
            }
        }
    }
    return output;
}

//This function searches for the best chi-value in a 3D besT_PSFs vector, and returns the location
vector<int> FindBestPSF(vector<vector<vector<double>>> best_PSFs){
    vector<int> best_fit{0,0,0};
    for (int i=0; i< best_PSFs.size();i++){
        for (int j=0; j< best_PSFs[0].size();j++){ 
            for (int k=0; k< best_PSFs[0][0].size();k++){
                double curr_best = best_PSFs[best_fit[0]][best_fit[1]][best_fit[2]];
                double curr_val = best_PSFs[i][j][k];
                if (curr_val < curr_best){
                    best_fit[0] = i;
                    best_fit[1] = j;
                    best_fit[2] = k;
                }
            }
        }       
    }
    
    return best_fit;
}
/*
double* Reshape2D(vector<double> chi_array, int shape1, int shape2){
    vector<vector<double>> output(shape1,shape2);
    for(i=0;i<shape1;i++){
        for(j=0;j<shape2;j++){
            output[i][j] = chi_array[i+j]
        }
    }
    return output;
}
*/

//This function Multiplies every value in a 2D double vector by a given double value 
vector<vector<double>> MultiplyByScalar2D(vector<vector<double>> vec, double val){
    vector<vector<double>> output(vec.size());
    for (int i=0;i<vec.size();i++){
        vector<double> row(vec[0].size());
        for (int j=0;j<vec[0].size();j++){
            row[j] =(vec[i][j]*val);
            //output[i][j] = vec[i][j]*val;
        }
        output[i] = row;
    }
    return output;
}

//This function Multiplies every value in a 1D double vector by a given double value 
vector<double> MultiplyByScalar1D(vector<double> vec, double val){
    vector<double> output(vec.size());
    for (int i=0;i<vec.size();i++){
        output[i] = vec[i]*val;
    }
    return output;
}

//This function Multiplies every value in a 1D int vector by a given double value 
vector<int> MultiplyByScalar1DInt(vector<int> vec, double val){
    vector<int> output(vec.size());
    for (int i=0;i<vec.size();i++){
        output[i] = (vec[i]*val);
    }
    return output;
}

//This function Sums the Square of all values in a 2D vector
double SumSquared2D(vector<vector<double>> vec1){
    double output = 0;
    for(int i =0; i< vec1.size();i++){
        for(int j=0; j< vec1[0].size();j++){
            output = output+(vec1[i][j]*vec1[i][j]);
        }
    }
    return output;
}

//This function Sums the Square of all values in a 1D double vector
double SumSquared1D(vector<double> vec1){
    double output = 0;
    for(int i =0; i< vec1.size();i++){
        output = output+(vec1[i]*vec1[i]);
    }
    return output;
}

//This function Sums the Square of all values in a 1D int vector
int SumSquared1DInt(vector<int> vec1){
    int output = 0;
    for(int i =0; i< vec1.size();i++){
        output = output+(vec1[i]*vec1[i]);
    }
    return output;
}

//This function Sums all of the values in a double 1D vector
double SumVector(vector<double> vec1){
    double output = 0;
    for(int i =0; i< vec1.size();i++){
        output = output+vec1[i];
    }
    return output;
}

//This function Sums all of the values in a int 1D vector
int SumVectorInt(vector<int> vec1){
    int output = 0;
    for(int i =0; i< vec1.size();i++){
        output = output+vec1[i];
    }
    return output;
}

//This function Adds all of the values in 2 2D vectors of the same size to each other, outputting a 2D double vector
vector<vector<double>> Add2DVectors(vector<vector<double>> vec1, vector<vector<double>> vec2){
    vector<vector<double>> output(vec1.size(), vector<double>(vec1[0].size(),0));
    for(int i=0;i<vec1.size();i++){
        vector<double> row(vec1[0].size());
        for (int j=0;j<vec1[0].size();j++){
            row[j] = (vec1[i][j]+vec2[i][j]);
        }
        output[i] = row;
    }
    return output;
}

//This function Adds all of the values in 2 1D vectors of the same size to each other, outputting a 1D double vector
vector<double> Add1DVectors(vector<double> vec1, vector<double> vec2){
    vector<double> output(vec1.size());
    for(int i=0;i<vec1.size();i++){
        output[i] = (vec1[i]+vec2[i]);
    }
    return output;
}

//This function multiplies all of the values in 2 2D vectors of the same size to each other, outputting a 2D double vector
vector<vector<double>> Multiply2DVectors(vector<vector<double>> vec1, vector<vector<double>> vec2){
    vector<vector<double>> output(vec1.size(), vector<double>(vec1[0].size(),0));
    for(int i=0;i<vec1.size();i++){
        vector<double> row(vec1[0].size());
        for (int j=0;j<vec1[0].size();j++){
            row[j] = (vec1[i][j]*vec2[i][j]);
        }
        output[i] = row;
    }
    return output;
}

//This function multiplies all of the values in 2 1D vectors (double,double) of the same size to each other, outputting a 1D double vector
vector<double> Multiply1DVectors(vector<double> vec1, vector<double> vec2){
    vector<double> output(vec1.size());
    for(int i=0;i<vec1.size();i++){       
        output[i] = vec1[i]*vec2[i];
    }
    return output;
}

//This function multiplies all of the values in 2 1D vectors (int, double) of the same size to each other, outputting a 1D double vector
vector<double> Multiply1DVectorsInt(vector<int> vec1, vector<double> vec2){
    vector<double> output(vec1.size());
    for(int i=0;i<vec1.size();i++){       
        output[i] = (double)vec1[i]*vec2[i];
    }
    return output;
}

//This function multiplies all of the values in 2 1D vectors (int, int) of the same size to each other, outputting a 1D double vector
vector<double> Multiply1DVectorsIntInt(vector<int> vec1, vector<int> vec2){
    vector<double> output(vec1.size());
    for(int i=0;i<vec1.size();i++){       
        output[i] = (double)vec1[i]*(double)vec2[i];
    }
    return output;
}

//This function subtracts all of the values in 2 2D vectors (double, double) of the same size to each other, outputting a 2D double vector
vector<vector<double>> Subtract2DVectors(vector<vector<double>> vec1, vector<vector<double>> vec2){
    vector<vector<double>> output(vec1.size(),vector<double>(vec1[0].size(),0));
    for(int i=0;i<vec1.size();i++){
        vector<double> row(vec1[0].size());
        for (int j=0;j<vec1[0].size();j++){
            row[j] = vec1[i][j]-vec2[i][j];
        }
        output[i] = row;
    }
    return output;
}

//This function subtracts all of the values in 2 1D vectors (double, double) of the same size to each other, outputting a 1D double vector
vector<double> Subtract1DVectors(vector<double> vec1, vector<double> vec2){
    vector<double> output(vec1.size());
    for(int i=0;i<vec1.size();i++){
        output[i] = vec1[i]-vec2[i];
    }
    return output;
}

//This function subtracts all of the values in 2 1D vectors (double, int) of the same size to each other, outputting a 2D double vector
vector<double> Subtract1DVectorsDoubInt(vector<double> vec1, vector<int> vec2){
    vector<double> output(vec1.size());
    for(int i=0;i<vec1.size();i++){
        output[i] = vec1[i]-(double)vec2[i];
    }
    return output;
}

//This function creates a binary PSF vecot rof ints describing where in the array of psfs the best fit psf is located.
vector<int> BinaryPSF(vector<vector<double>> array_of_PSFs, pair<int,int> secondary_coordinates, vector<double> image_psf, double flux, vector<double> w)
{
    //cout << secondary_coordinates << endl;
    
	vector<double> scale_array;
    scale_array.reserve(19);
    for(double l= 0.95; l > 0; l = l-0.05){
        scale_array.push_back(l);
    }
    
	int models_length = array_of_PSFs.size();
    cout << "modles length " << models_length << endl;
	int scale_length = scale_array.size();
	int chi_array_length = models_length * models_length * scale_length;
	
	vector<double> chi_array(chi_array_length);
	
	int cx = (array_of_PSFs.size() - 1) * .5;
	int cy = (array_of_PSFs[0].size() - 1) * .5;

	//cout << "line 303" << endl;
	for (int i = 0; i < scale_length; i++) { 
        double relative_flux = scale_array[i];
        //cout << "line 332" << endl;
		
        //A slice of model_PSF from model_PSF[cy-2:cy+3,cx-2:cx+3]
		//vector<double> PSF1 = (xt::view(model_PSF, xt::range(cy - 2, cy + 3 ), xt::range(cx - 2, cx + 3)));
        //cout << "361 " << endl;
		for (int j = 0; j < models_length; j++) {
            //cout << "line 348" << endl;
			vector<double> model_PSF = array_of_PSFs[j];            
            int rangeY_1 = cy - secondary_coordinates.first;
			int rangeY_2 = cy - secondary_coordinates.first + 5;
			int rangeX_1 = cx - secondary_coordinates.second;
			int rangeX_2 = cx - secondary_coordinates.second + 5;
            vector<double> PSF2((rangeY_2-rangeY_1)+(rangeX_2-rangeX_1));
            int n = 0;
            for (int m = rangeY_1; m < rangeY_2; m++){
                PSF2[n] = model_PSF[m];
                n++;
            }
            for (int m = rangeX_1; m < rangeX_2; m++){
                PSF2[n] = model_PSF[m];
                n++;
            }
            
            //vector<double> a = PSF1;
            //vector<double> b = PSF2;
            //cout << "377 " << scale_length << endl;
			for (int k = 0; k < models_length; k++) {    
                vector<double> model_PSF = array_of_PSFs[k];
                vector<double> PSF1(10);
                n = 0;
                for (int b = cy-2; b < cy+3; b++){
                    PSF1[n] = model_PSF[b];
                    n++;
                }
                for (int b = cx-2; b < cx+3; b++){
                    PSF1[n] = model_PSF[b];
                    n++;
                }
				
                
                //a = MultiplyByScalar1D(a, relative_flux);
                //b = MultiplyByScalar1D(b, (1-relative_flux));
                //vector<double> PSF_sum = Add1DVectors(MultiplyByScalar1D(PSF1, relative_flux), MultiplyByScalar1D(PSF2, (1-relative_flux)));
                //vector<double> temp_array = MultiplyByScalar1D(PSF_sum,flux);
                //cout << "line 375" << endl;
                vector<double> chi_PSF = Subtract1DVectors(image_psf,MultiplyByScalar1D(Add1DVectors(MultiplyByScalar1D(PSF1, relative_flux), MultiplyByScalar1D(PSF2, (1-relative_flux))),flux));
                //chi_PSF = Multiply1DVectors(chi_PSF, w);
				//vector<double> chi_PSF = w * (image_psf - temp_array);
                //cout << "line 379" << endl;
                double chi_value = SumSquared1D(Multiply1DVectors(chi_PSF, w));
				//double chi_value = sum(chi_PSF * chi_PSF)(0);
                //cout << "line 382 " << i*models_length*models_length+j*models_length+k << endl;
				chi_array[i*models_length*models_length+j*models_length+k] = chi_value;	
                //cout << "line 385 " << k << endl;
			}
            //cout << "line 387 " << j << endl;
		}
        //cout << "line 389 " << i << endl;
	}
    cout << "Hi" << endl;
    cout << "Hey there " << chi_array.size() << endl;
    vector<vector<vector<double>>> chi = Reshape3D(chi_array, 100,100,19);

    vector<int> best_PSFs{0,0,0};
    //This loops searches for the best chi value, and saves the best_PSFs as the indices that equate the best fit
    for (int i=0; i<chi.size();i++){
        for (int j=0; j<chi[i].size();j++){
            for (int k=0; k<chi[i][j].size();k++){
                double new_val = chi[i][j][k];
                
                double best = chi[best_PSFs[0]][best_PSFs[1]][best_PSFs[2]];
                if (new_val < best) { 
                    best_PSFs[0] = i;
                    best_PSFs[1] = j;
                    best_PSFs[2] = k;
                }
            }
        }
    }

	return best_PSFs;
}

//This function returns a vector of the best relative fluxs for the best binary fit given
vector<double> BinaryRelativeFlux(vector<vector<double>> array_of_PSFs, pair<int,int> secondary_coordinates, vector<double> image_psf, double flux, vector<int> w, int best_primary, int best_secondary) 
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
	
	primary_ref[0] = best_primary - 11;
	primary_ref[1] = best_primary - 10;
	primary_ref[2] = best_primary - 9;
	primary_ref[3] = best_primary - 1;
	primary_ref[4] = best_primary;
	primary_ref[5] = best_primary + 1;
	primary_ref[6] = best_primary + 9;
	primary_ref[7] = best_primary + 10;
	primary_ref[8] = best_primary + 11;
	
	secondary_ref[0] = best_secondary - 11;
	secondary_ref[1] = best_secondary - 10;
	secondary_ref[2] = best_secondary - 9;
	secondary_ref[3] = best_secondary - 1;
	secondary_ref[4] = best_secondary;
	secondary_ref[5] = best_secondary + 1;
	secondary_ref[6] = best_secondary + 9;
	secondary_ref[7] = best_secondary + 10;
	secondary_ref[8] = best_secondary + 11;
	vector<double> scale_array;
    scale_array.reserve(99);
    for(double i = 0.99; i > 0; i=i-0.01){
        scale_array.push_back(i);
    }
	//vector<double> scale_array = xt::arange<double>(.99, 0, -.01);
	
	int primary_length = primary_ref.size();
	int secondary_length = secondary_ref.size();
	int scale_length = scale_array.size();
	int chi_array_length = primary_length * secondary_length * scale_length;
	vector<double> chi_array(chi_array_length);
	
	for (int i = 0; i < primary_length; ++i ) {
		int primary = primary_ref[i];
		vector<double> model_PSF = array_of_PSFs[primary];
		int cx = (model_PSF.size() - 1) / 2;
		int cy = (model_PSF.size() - 1) / 2;

        vector<double> PSF1(10);
        int k =0 ;
        for (int j = (cy-2); j < (cy+3); j++){
            PSF1[j] = model_PSF[j];
            k++;
        }
        for (int j = (cx-2); j < (cx+3); j++){
            PSF1[j] = model_PSF[j];
            k++;
        }
        //PSF1.push_back(cy_range);
        //PSF1.push_back(cx_range);
		//vector<double> PSF1 = (xt::view(model_PSF, xt::range(cy - 2, cy + 3 ), xt::range(cx - 2, cx + 3)));	
		for (int j = 0; j < secondary_length; ++j) {
			int secondary = secondary_ref[j];
			vector<double> model_PSF = array_of_PSFs[secondary];
			int rangeY_1 = cy - secondary_coordinates.first;
			int rangeY_2 = cy - secondary_coordinates.first + 5;
			int rangeX_1 = cx - secondary_coordinates.second;
			int rangeX_2 = cx - secondary_coordinates.second + 5;
			

            vector<double> PSF2((rangeY_2-rangeY_1)+(rangeX_2-rangeX_1));
            int k = 0;
            for (int n = rangeY_1; n < rangeY_2; n++){
                PSF2[k] = model_PSF[n];
                k++;
            }
            for (int n = rangeX_1; n < rangeX_2; n++){
                PSF2[k] = model_PSF[n];
                k++;
            }

			//vector<double> PSF2 = (xt::view(model_PSF, xt::range(rangeY_1, rangeY_2), xt::range(rangeX_1, rangeX_2)));	
			for (int k = 0; k < scale_length; ++k) {
				double relative_flux = scale_array[k];
				vector<double> a = MultiplyByScalar1D(PSF1, relative_flux);
				vector<double> b = MultiplyByScalar1D(PSF2, (1 - relative_flux));
				vector<double> temp_array = MultiplyByScalar1D(Add1DVectors(a, b), flux);
                
                vector<double> temp_sub = Subtract1DVectors(image_psf , temp_array);
                
				vector<double> chi_PSF = Multiply1DVectorsInt(w, temp_sub);
                
				double chi_value = SumSquared1D(chi_PSF);
				chi_array[i*primary_length*secondary_length+j*primary_length+k] = chi_value;
			}
		}
	}
	
	//vector<double> chi = xt::adapt(chi_array, {9, 9, 99});
    //for (double val: chi_array) { cout << val << endl; } 
    vector<vector<vector<double>>> chi = Reshape3D(chi_array, 9, 9, 99);
    vector<int> best_PSFs = FindBestPSF(chi);
	//vector<int> best_PSFs = xt::adapt(xt::unravel_index(xt::argmin(chi)(0), chi.size(), xt::layout_type::row_major), {3});
	double best_primary_ref = primary_ref[best_PSFs[0]];
	double best_secondary_ref = secondary_ref[best_PSFs[1]];
	double best_ref = scale_array[best_PSFs[2]];
	vector<double> final_best_PSF = {best_primary_ref, best_secondary_ref, best_ref};
	
	return final_best_PSF;
}

//This function returns the best flux value of the BinaryFit for the run output
double BinaryFit(vector<double> model_primary, vector<double> model_secondary, pair<int,int> secondary_coordinates, vector<double> image_psf, double base, vector<int> w, double relative_flux) 
{
	int cx = (model_primary.size() - 1) / 2;
	int cy = (model_primary.size() - 1) / 2;
	int rangeY_1 = cy - secondary_coordinates.first;
	int rangeY_2 = cy - secondary_coordinates.first + 5;
	int rangeX_1 = cx - secondary_coordinates.second;
	int rangeX_2 = cx - secondary_coordinates.second + 5;
	
    vector<double> PSF1(10);
    int k=0;
    for (int j = (cy-2); j < (cy+3); j++){
        //cout << j;
        PSF1[k] = model_primary[j];
        k++;
    }
    for (int j = (cx-2); j < (cx+3); j++){
        //cout << j;
        PSF1[k] = model_primary[j];
        k++;
    }
    //PSF1.push_back(cy_range);
    //PSF1.push_back(cx_range);
    vector<double> PSF2((rangeY_2-rangeY_1)+(rangeX_2-rangeX_1));
    k = 0;
    for (int j = (rangeY_1); j < (rangeY_2); j++){
        //cout << j;
        PSF2[k] = model_secondary[j];
        k++;
    }
    for (int j = (rangeX_1); j < (rangeX_2); j++){
        //cout << j;
        PSF2[k] = model_secondary[j];
        k++;
    }
	vector<double> flux_array;
	vector<double> chi_array(1001);
	
	//change "range" to something like "Scale"
    
	vector<double> scale_array;
    scale_array.reserve(1001);
    // = xt::arange<double>(0, 10.01, 0.01);
    int j=0;
    for (double i=0.0;i<10.01;i=i+0.01){
        scale_array.push_back(i);
    }
	int scale_length = scale_array.size();
	flux_array.reserve(scale_length);
	chi_array.reserve(scale_length);
    vector<double> a = MultiplyByScalar1D(PSF1, relative_flux);
    vector<double> b =  MultiplyByScalar1D(PSF2, (1 - relative_flux));
	vector<double> PSF_sum = Add1DVectors(a,b);

    
	for (int i = 0; i < scale_length; ++i) {
		double scale = scale_array[i];
        flux_array.push_back(scale * base);
        
		vector<double> scale_PSF = MultiplyByScalar1D(PSF_sum, scale*base);
        //cout << scale_PSF.size() << endl;
        
        
		vector<double> chi_PSF = Multiply1DVectorsInt(w, Subtract1DVectors(image_psf, scale_PSF));
		double chi_value = SumSquared1D(chi_PSF);
		chi_array[i] = chi_value;
	}
		
	//auto chi_squareds = xt::adapt(chi_array, {scale_length});
	//auto fluxes = xt::adapt(flux_array, {scale_length});
    
	double min_chi = chi_array[0];
    for(int i=1;i<chi_array.size();i++){
        if (chi_array[i] < min_chi) { min_chi = chi_array[i]; }
    }
	int min_chi_index = 0;
	int chi_length = chi_array.size();
	//USE A VARIABLE FOR CHI_SQUAREDS.SHAPE, CHECK IN THE OTHER FILES TOO
	for (int i = 0; i < chi_length; ++i) {
		if (chi_array[i] == min_chi) {
			min_chi_index = i;
			break;
		}
	}	
	
	double best_flux = flux_array[min_chi_index];
	
	return best_flux;	 
}

//This function returns the flux number of the secondary object in the binary fit
double SecondaryBinaryFit(vector<double> primary, vector<double> secondary, pair<int,int> secondary_coordinates, vector<int> w2, vector<double> image_psf, double flux, double relative_flux)
{
	int cx = (primary.size() - 1) / 2;
	int cy = (primary.size() - 1) / 2;
	int rangeY_1 = cy - secondary_coordinates.first;
	int rangeY_2 = cy - secondary_coordinates.first + 5;
	int rangeX_1 = cx - secondary_coordinates.second;
	int rangeX_2 = cx - secondary_coordinates.second + 5;
	vector<double> PSF1(10);
    vector<double> PSF2((rangeY_2-rangeY_1)+(rangeX_1-rangeX_1));
    int k=0;
    for (int j = cy-2; j < cy+3; j++){
        PSF1[k] = primary[j];
        k++;
    }
    for (int j = cx-2; j < cx+3; j++){
        PSF1[k] = primary[j];
        k++;
    }
    //PSF1.push_back(cy_range);
    //PSF1.push_back(cx_range);
    k = 0;
    for (int j = rangeY_1; j < rangeY_2; j++){
        PSF2[k] = secondary[j];
        k++;
    }
    for (int j = rangeX_1; j < rangeX_2; j++){
        PSF2[k] = secondary[j];
        k++;
    }
    
    //cout << PSF1.size() << endl;
    //cout << PSF2.size() << endl;
	vector<double> scale_array(100001);
    //int j = 0;
    for(double i=0;i<100001; i++){
        scale_array.push_back(i*0.0001);
    }
    
    
	int scale_length = scale_array.size();
	vector<double> flux_array;
	vector<double> chi_array(scale_length);
	flux_array.reserve(scale_length);
	chi_array.reserve(scale_length);
    
	vector<double> temp1 = MultiplyByScalar1D(PSF1, flux*relative_flux);
    vector<double> temp2 = Subtract1DVectors(image_psf, temp1);
	vector<double> relative_PSF2 = Multiply1DVectorsInt(w2, temp2);
                                                     
                                                     
	for (int i = 0; i < scale_length; ++i)
	{
		double scale = scale_array[i];
		vector<double> scale_PSF2 =  MultiplyByScalar1D(MultiplyByScalar1D(Multiply1DVectorsInt(w2,PSF2),scale),flux);
		flux_array.push_back(flux* scale);
		vector<double> chi_PSF = Multiply1DVectorsInt(w2,Subtract1DVectors(relative_PSF2,scale_PSF2));
		double chi_value = SumSquared1D(chi_PSF);
		chi_array[i] = chi_value;
	}

	//auto chi_squareds = xt::adapt(chi_array, {scale_length});
	//auto fluxes = xt::adapt(flux_array, {scale_length});
    double min_chi = chi_array[0];
    for(int i=1;i<chi_array.size();i++){
        if (chi_array[i] < min_chi) { min_chi = chi_array[i]; }
    }
	//double min_chi = xt::amin(chi_squareds)(0);
	int min_chi_index = 0;
	int chi_length = chi_array.size();
	for (int i = 0; i < chi_length ; ++i) {
		if (chi_array[i] == min_chi) {
			min_chi_index = i;
			break;
		}
	}	                                                    
	
	double best_flux = flux_array[min_chi_index];
	
	return best_flux;
}

                                                    
                                                     
//This function finds the relatvie ebst flux for the primary and secondary objects in the fit
double RelativeFluxes(vector<double> PSF_1, pair<int,int> coordinates_1, vector<double> PSF_2, pair<int,int> coordinates_2, vector<int> w, vector<double> image_psf, double base, int flux) 
{
	int cx = (PSF_1.size() - 1) * .5;
	int cy = (PSF_2.size() - 1) * .5;
	
	int P_rangeY_1 = cy - coordinates_1.first;
	int P_rangeY_2 = cy - coordinates_1.first + 5;
	int P_rangeX_1 = cx - coordinates_1.second;
	int P_rangeX_2 = cx - coordinates_1.second + 5;
	
	int S_rangeY_1 = cy - coordinates_2.first;
	int S_rangeY_2 = cy - coordinates_2.first + 5;
	int S_rangeX_1 = cx - coordinates_2.second;
	int S_rangeX_2 = cx - coordinates_2.second + 5;

    vector<int> PSF1((P_rangeY_2-P_rangeY_1)+(P_rangeX_2-P_rangeX_1));
    int k=0;
    for (int j = P_rangeY_1; j < P_rangeY_2; j++){
        PSF1[k] = PSF_1[j];
        k++;
    }
    for (int j = P_rangeX_1; j < P_rangeX_2; j++){
        PSF1[k] = PSF_1[j];
        k++;
    }
    //PSF1.push_back(cy_range);
    //PSF1.push_back(cx_range);
    vector<int> PSF2((S_rangeY_2-S_rangeY_1)+(S_rangeX_2-S_rangeX_1));
    k=0;
    for (int j = S_rangeY_1; j < S_rangeY_2; j++){
        PSF2[k] = PSF_2[j];
        k++;
    }
    for (int j = S_rangeX_1; j < S_rangeX_2; j++){
        PSF2[k] = PSF_2[j];
        k++;
    }
    //PSF2.push_back(cy_range);
    //PSF2.push_back(cx_range);
	//vector<double> PSF1 = (xt::view(PSF_1, xt::range(P_rangeY_1, P_rangeY_2), xt::range(P_rangeX_1, P_rangeX_2)));
	//vector<double> PSF2 = (xt::view(PSF_2, xt::range(S_rangeY_1, S_rangeY_2), xt::range(S_rangeX_1, S_rangeX_2)));
	
	vector<double> scale_array(100001);// = xt::arange<double>(0, 10.0001, .0001);
    int j = 0;
    for(double i=0;i<100001; i++){
        scale_array.push_back(i*0.0001);
    }
	int scale_length = scale_array.size();
	vector<double> flux_array;
	vector<double> chi_array(scale_length);
	flux_array.reserve(scale_length);
	chi_array.reserve(scale_length);
	
	vector<double> rPSF = Multiply1DVectorsInt(w, Subtract1DVectorsDoubInt(image_psf, MultiplyByScalar1DInt(PSF2, flux)));
	for (int i = 0; i < scale_length; ++i) {
		double scale = scale_array[i];
		vector<double> scale_PSF = MultiplyByScalar1D(MultiplyByScalar1D(Multiply1DVectorsIntInt(PSF1,w),scale),base);
		flux_array.push_back(scale * base);
		vector<double> chi_PSF = Multiply1DVectorsInt(w, Subtract1DVectors(rPSF,scale_PSF));
		double chi_value = SumSquared1D(chi_PSF);
		chi_array[i] = chi_value;
	}
	
	//auto chi_squareds = xt::adapt(chi_array, {scale_length});
	//auto fluxes = xt::adapt(flux_array, {scale_length});
	//double min_chi = xt::amin(chi_squareds)(0);
    double min_chi = chi_array[0];
    for(int i=1;i<chi_array.size();i++){
        if (chi_array[i] < min_chi) { min_chi = chi_array[i]; }
    }
	int min_chi_index = 0;
	int chi_length = chi_array.size();
	for (int i = 0; i < chi_length ; ++i) {
		if (chi_array[i] == min_chi) {
			min_chi_index = i;
			break;
		}
	}	
	
	double best_flux = flux_array[min_chi_index];
	
	return best_flux;
}

//This function reads in an int txt file of a vecotr of pairs
vector<pair<int,int>> Read_FileInt(string filename){
    vector<pair<int,int>> vec; 
    short loop=0; 
    //vector<double> line; 
    ifstream myfile(filename); 
    double val1, val2;
    while (!myfile.eof()){
        myfile >> val1 >> val2;
        
        pair<int,int> temp_pair((int)val1,(int)val2);
        //cout << temp_pair.first << " " << temp_pair.second << endl;
        vec.push_back(temp_pair); 
    }
    
    vec.pop_back();

    return vec;
}
        
//This function reads in an double txt file of a vecotr of pairs
vector<pair<int,int>> Read_FileDouble(string filename){
    //vector<vector<double>> vec; 
    vector<pair<int,int>> vec;
    short loop=0; 
    //vector<double> line; 
    ifstream myfile(filename); 
    double val1, val2;
    if (myfile.is_open()){
        while (!myfile.eof()){
            myfile >> val1 >> val2;
        
            pair<int,int> temp_pair((int)val1,(int)val2);
        
            vec.push_back(temp_pair); 
        }
        myfile.close();
    }
    return vec;
}

//This function reads a txt file into a vecctor of doubles
vector<double> Read_FileDouble1D(string filename){
    vector<double> vec; 

    ifstream myfile(filename); 
    int row;
    if (myfile.is_open()) 
    {
        while(!myfile.eof()){
            double val1;
            myfile >> val1;
            vec.push_back(val1);
        }  
        myfile.close();
    }
    //for(double val : vec){ cout << val << " ";}
    //cout << endl;
    return vec;
}



//MAIN 
int main() {    
    //read in the image and model data, and the parameters we set earlier
    string ipsf = "tempdata/image_psf.txt";
    vector<double> image_psf = Read_FileDouble1D(ipsf); 
    //for(auto val : image_psf){ cout << val << " ";}

    int array_length = 100;
    vector<vector<double>> array_of_PSFs(100);
    
    for (int i = 0; i <= array_length; ++i) {
    		string filename = "tempdata/PSFmodel_" + to_string(i) + ".txt";
    		array_of_PSFs[i - 1] = Read_FileDouble1D(filename);
    	}
    
    	
    string seccoords = "tempdata/sec_coords.txt";
    string primcoords = "tempdata/prim_coords.txt";
	vector<pair<int,int>> coordinates = Read_FileInt(seccoords);
	vector<pair<int,int>> center_load = Read_FileInt(primcoords);
	pair<int,int> center = center_load[0];
    	
    	pair<int,int> primary_coordinates{2,2};
    	pair<int,int> secondary_coordinates;
    	
    	//arrays to store the output values
    	vector<double> residual_error_array(coordinates.size());
	vector<int> best_primary_array(coordinates.size());
	vector<double> flux_primary_array(coordinates.size());
	vector<int> best_secondary_array(coordinates.size());
	vector<double> flux_secondary_array(coordinates.size());
	vector<double> flux_sum_array(coordinates.size());

        vector<pair<int,int>> temp_coords{make_pair(1,1),make_pair(1,2),make_pair(1,3),make_pair(2,1),make_pair(2,2),make_pair(2,3),make_pair(3,1),make_pair(3,2),make_pair(3,3)};
    
        //This for loop tests every pixel of the 9 surrounding the centroid
    for (int i = 0; i < coordinates.size(); i++) {

    		pair<int,int> secondary = coordinates[i];

            secondary_coordinates = temp_coords[i];
            
    		vector<vector<int>> w_array = BinaryBase(primary_coordinates, secondary_coordinates);

            vector<int> w = w_array[0];       
            vector<int> w1 = w_array[1];       
            vector<int> w2 = w_array[2];       
            
            vector<double> sqrt_weights;
            vector<int> inted_sqrtw;
            for(int j=0;j<image_psf.size();j++){
                double val = image_psf[j];
                if (val<0){ val = val*-1; }
                sqrt_weights.push_back(sqrt(val));
                inted_sqrtw.push_back((int)sqrt(val));
            }

    		int iteration = 0;
    		int best_primary = 0;
    		int best_secondary = 0;
    		double best_RF = 0.0;
            
    		double base = SumVector(Multiply1DVectorsInt(w, image_psf));

    		double comp = 0;
    		double flux = base * 1.0;

    		while (abs(comp - flux) > .000001 && iteration < 10) {
    			comp = flux;
                //cout << "line 398" << endl;
                //Find the Binary PSF indices
    			vector<int> best_PSFs =  BinaryPSF(array_of_PSFs, secondary_coordinates, image_psf, flux, sqrt_weights);
                //cout << "line 401" << endl;
    			best_primary = best_PSFs[0];
    			best_secondary = best_PSFs[1];
                
                //Find the relative flux related to the indices
    			vector<double> best_PSFs_ref = BinaryRelativeFlux(array_of_PSFs, secondary_coordinates, image_psf, flux, inted_sqrtw, best_primary, best_secondary);
    			best_primary = best_PSFs_ref[0];
    			best_secondary = best_PSFs_ref[1];
    			best_RF = best_PSFs_ref[2];
    			//cout << "line 924" << endl;
			flux = BinaryFit(array_of_PSFs[best_primary], array_of_PSFs[best_secondary], secondary_coordinates, image_psf, base, w, best_RF);

    			++iteration;
    		}
            
    		

    		vector<double> primary_PSF = array_of_PSFs[best_primary];
            
    		vector<double> secondary_PSF = array_of_PSFs[best_secondary];

    		double comp_primary = 1.0;
    		double comp_secondary = 1.0;
    		double flux_primary = 2.0;
            //cout << "line 936" << endl;
    		double flux_secondary = SecondaryBinaryFit(primary_PSF, secondary_PSF, secondary_coordinates, w2, image_psf, flux, best_RF);
            

    		iteration = 0;
            
            
            while (abs(comp_primary - flux_primary) > .000001 && abs(comp_secondary - flux_secondary) > .000001) {
                comp_primary = flux_primary;
                comp_secondary = flux_secondary;
                if (iteration >= 10) {break;}
                flux_primary = RelativeFluxes(primary_PSF, primary_coordinates, secondary_PSF, secondary_coordinates, w1, image_psf, base, flux_secondary);
                flux_secondary = RelativeFluxes(secondary_PSF, secondary_coordinates, primary_PSF, primary_coordinates, w2, image_psf, base, flux_primary);

                ++iteration;
            }


		int cx = (primary_PSF.size() - 1) * .5;
		int cy = (primary_PSF.size() - 1) * .5;
		
		int P_rangeY_1 = cy - primary_coordinates.first;
		int P_rangeY_2 = cy - primary_coordinates.first + 5;
		int P_rangeX_1 = cx - primary_coordinates.second;
		int P_rangeX_2 = cx - primary_coordinates.second + 5;
		
		int S_rangeY_1 = cy - secondary_coordinates.first;
		int S_rangeY_2 = cy - secondary_coordinates.first + 5;
		int S_rangeX_1 = cx - secondary_coordinates.second;
		int S_rangeX_2 = cx - secondary_coordinates.second + 5;
            
       
        vector<double> PSF1;
        vector<double> PSF2;
        
        for (int j = P_rangeY_1; j < P_rangeY_2; j++) {
            PSF1.push_back(primary_PSF[j]);
        }
        for (int j = S_rangeY_1; j < S_rangeY_2; j++) {
            PSF2.push_back(primary_PSF[j]);
        }
		vector<double> PSF_sum = Add1DVectors(MultiplyByScalar1D(PSF1,flux_primary), MultiplyByScalar1D(PSF2,flux_secondary));
            
		vector<double> pre_residual = Subtract1DVectors(image_psf,PSF_sum);
        for(int j=0;j< pre_residual.size();j++){
            if(pre_residual[j] < 0) { pre_residual[j] = pre_residual[j]*-1; }
        }
        vector<double> temp_vec = Multiply1DVectorsInt(w,pre_residual);
		double residual = SumVector(temp_vec);
		double residual_error = residual / SumVectorInt(w);
		
		cout << "psf1: " << best_primary << " flux1: " << flux_primary << " psf2: " << best_secondary << " flux2: " << flux_secondary << " error: " << residual_error << endl;

		residual_error_array[i] = residual_error;
		best_primary_array[i] = best_primary;
		flux_primary_array[i] = flux_primary;
		best_secondary_array[i] = best_secondary;
		flux_secondary_array[i] = flux_secondary;
		flux_sum_array[i] = flux_primary + flux_secondary;
    }
            
    
	
	cout << endl;
    
            
    pair<int,int> cent(center.first,center.second);
    	vector<pair<int,int>> center_array{cent, cent, cent, cent, cent, cent, cent, cent, cent};
    	vector<pair<int,int>> secondary_array = coordinates;
    	
    	//send the "outputs" as seperate files since they are different data types and we can read them all in to python for the calculation of angles and separation - it only takes a second :)
    ofstream outFile1("tempdata/residual_error_array.txt");
    for (const auto &e : residual_error_array) outFile1 << e << "\n";
    ofstream outFile2("tempdata/center_array.txt");
    for (const auto &e : center_array){ 
        outFile2 << e.first << " " << e.second << "\n";
    }
    ofstream outFile3("tempdata/secondary_array.txt");
    for (const auto &e : secondary_array){ 
            outFile3 << e.first << " " << e.second << "\n";
    }
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
