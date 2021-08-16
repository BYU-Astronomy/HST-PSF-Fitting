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

vector<vector<int>> BinaryBase(vector<int> primary_coordinates, vector<int> secondary_coordinates) 
{
	vector<vector<int>> w1 {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
	vector<vector<int>> w2 {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
    vector<vector<int>> w{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};   
	w1[primary_coordinates[0] - 1][ primary_coordinates[1] - 1] = 1;
	w1[primary_coordinates[0] - 1][ primary_coordinates[1]] = 1;
	w1[primary_coordinates[0] - 1][ primary_coordinates[1] + 1] = 1;
	w1[primary_coordinates[0]][ primary_coordinates[1] - 1] = 1;
	w1[primary_coordinates[0]][ primary_coordinates[1]] = 1;
	w1[primary_coordinates[0]][ primary_coordinates[1] + 1] = 1;
	w1[primary_coordinates[0] + 1][ primary_coordinates[1] - 1] = 1;
	w1[primary_coordinates[0] + 1][ primary_coordinates[1]] = 1;
	w1[primary_coordinates[0] + 1][ primary_coordinates[1] + 1] = 1;

    w[primary_coordinates[0] - 1][ primary_coordinates[1] - 1] = 1;
	w[primary_coordinates[0] - 1][ primary_coordinates[1]] = 1;
	w[primary_coordinates[0] - 1][ primary_coordinates[1] + 1] = 1;
	w[primary_coordinates[0]][ primary_coordinates[1] - 1] = 1;
	w[primary_coordinates[0]][ primary_coordinates[1]] = 1;
	w[primary_coordinates[0]][ primary_coordinates[1] + 1] = 1;
	w[primary_coordinates[0] + 1][ primary_coordinates[1] - 1] = 1;
	w[primary_coordinates[0] + 1][ primary_coordinates[1]] = 1;
	w[primary_coordinates[0] + 1][ primary_coordinates[1] + 1] = 1;
    
	w2[secondary_coordinates[0] - 1][ secondary_coordinates[1] - 1] = 1;
	w2[secondary_coordinates[0] - 1][ secondary_coordinates[1]] = 1;
	w2[secondary_coordinates[0] - 1][ secondary_coordinates[1] + 1] = 1;
	w2[secondary_coordinates[0]][ secondary_coordinates[1] - 1] = 1;
	w2[secondary_coordinates[0]][ secondary_coordinates[1]] = 1;
	w2[secondary_coordinates[0]][ secondary_coordinates[1] + 1] = 1;
	w2[secondary_coordinates[0] + 1][ secondary_coordinates[1] - 1] = 1;
	w2[secondary_coordinates[0] + 1][ secondary_coordinates[1]] = 1;
	w2[secondary_coordinates[0] + 1][ secondary_coordinates[1] + 1] = 1;
    
	w[secondary_coordinates[0] - 1][ secondary_coordinates[1] - 1] = 1;
	w[secondary_coordinates[0] - 1][ secondary_coordinates[1]] = 1;
	w[secondary_coordinates[0] - 1][ secondary_coordinates[1] + 1] = 1;
	w[secondary_coordinates[0]][ secondary_coordinates[1] - 1] = 1;
	w[secondary_coordinates[0]][ secondary_coordinates[1]] = 1;
	w[secondary_coordinates[0]][ secondary_coordinates[1] + 1] = 1;
	w[secondary_coordinates[0] + 1][ secondary_coordinates[1] - 1] = 1;
	w[secondary_coordinates[0] + 1][ secondary_coordinates[1]] = 1;
	w[secondary_coordinates[0] + 1][ secondary_coordinates[1] + 1] = 1;
	

	vector<vector<int>> w_array{};
    w_array.insert( w_array.end(), w.begin(), w.end() );
    w_array.insert( w_array.end(), w1.begin(), w1.end() );
    w_array.insert( w_array.end(), w2.begin(), w2.end() );
	
	return w_array;
}

vector<vector<vector<double>>> Reshape3D(vector<double> chi_array, int shape1, int shape2, int shape3){
    vector<vector<vector<double>>> output(shape1,vector<vector<double>>(shape2,vector<double>(shape3)));
    for(int i=0;i<shape1;i++){
        for(int j=0;j<shape2;j++){
            for(int k=0;k<shape3;k++){
                output[i][j][k] = chi_array[i+j+k];
            }
        }
    }
    return output;
}

vector<vector<vector<int>>> Reshape3DInt(vector<double> chi_array, int shape1, int shape2, int shape3){
    vector<vector<vector<int>>> output(shape1,vector<vector<int>>(shape2,vector<int>(shape3)));
    for(int i=0;i<shape1;i++){
        for(int j=0;j<shape2;j++){
            for(int k=0;k<shape3;k++){
                output[i][j][k] = chi_array[i+j+k];
            }
        }
    }
    return output;
}

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
vector<vector<double>> MultiplyByScalar2D(vector<vector<double>> vec, double val){
    vector<vector<double>> output;
    for (int i=0;i<vec.size();i++){
        vector<double> row;
        for (int j=0;j<vec.size();j++){
            row.push_back(vec[i][j]*val);
            //output[i][j] = vec[i][j]*val;
        }
        output.push_back(row);
    }
    return output;
}

vector<double> MultiplyByScalar1D(vector<double> vec, double val){
    vector<double> output = vec;
    for (int i=0;i<vec.size();i++){
        output[i] = vec[i]*val;
    }
    return output;
}

vector<int> MultiplyByScalar1DInt(vector<int> vec, double val){
    vector<int> output = vec;
    for (int i=0;i<vec.size();i++){
        output[i] = vec[i]*val;
    }
    return output;
}

double SumSquared2D(vector<vector<double>> vec1){
    double output = 0;
    for(int i =0; i< vec1.size();i++){
        for(int j=0; j< vec1[0].size();j++){
            output = output+(vec1[i][j]*vec1[i][j]);
        }
    }
    return output;
}

double SumSquared1D(vector<double> vec1){
    double output = 0;
    for(int i =0; i< vec1.size();i++){
        output = output+(vec1[i]*vec1[i]);
    }
    return output;
}

int SumSquared1DInt(vector<int> vec1){
    int output = 0;
    for(int i =0; i< vec1.size();i++){
        output = output+(vec1[i]*vec1[i]);
    }
    return output;
}

double SumVector(vector<double> vec1){
    double output = 0;
    for(int i =0; i< vec1.size();i++){
        output = output+vec1[i];
    }
    return output;
}

int SumVectorInt(vector<int> vec1){
    int output = 0;
    for(int i =0; i< vec1.size();i++){
        output = output+vec1[i];
    }
    return output;
}

vector<vector<double>> Add2DVectors(vector<vector<double>> vec1, vector<vector<double>> vec2){
    vector<vector<double>> output;
    for(int i=0;i<vec1.size();i++){
        vector<double> row;
        for (int j=0;j<vec1[0].size();j++){
            row.push_back(vec1[i][j]+vec2[i][j]);
        }
        output.push_back(row);
    }
    return output;
}

vector<double> Add1DVectors(vector<double> vec1, vector<double> vec2){
    vector<double> output;
    for(int i=0;i<vec1.size();i++){
        output.push_back(vec1[i]+vec2[i]);
    }
    return output;
}

vector<vector<double>> Multiply2DVectors(vector<vector<double>> vec1, vector<vector<double>> vec2){
    vector<vector<double>> output;
    for(int i=0;i<vec1.size();i++){
        vector<double> row;
        for (int j=0;j<vec1[0].size();j++){
            row.push_back(vec1[i][j]*vec2[i][j]);
        }
        output.push_back(row);
    }
    return output;
}

vector<double> Multiply1DVectors(vector<double> vec1, vector<double> vec2){
    vector<double> output;
    for(int i=0;i<vec1.size();i++){       
        output.push_back(vec1[i]*vec2[i]);
    }
    return output;
}

vector<double> Multiply1DVectorsInt(vector<int> vec1, vector<double> vec2){
    vector<double> output;
    for(int i=0;i<vec1.size();i++){       
        output.push_back((double)vec1[i]*vec2[i]);
    }
    return output;
}

vector<double> Multiply1DVectorsIntInt(vector<int> vec1, vector<int> vec2){
    vector<double> output;
    for(int i=0;i<vec1.size();i++){       
        output.push_back((double)vec1[i]*(double)vec2[i]);
    }
    return output;
}

vector<vector<double>> Subtract2DVectors(vector<vector<double>> vec1, vector<vector<double>> vec2){
    vector<vector<double>> output;
    for(int i=0;i<vec1.size();i++){
        vector<double> row;
        for (int j=0;j<vec1[0].size();j++){
            row.push_back(vec1[i][j]-vec2[i][j]);
        }
        output.push_back(row);
    }
    return output;
}

vector<double> Subtract1DVectors(vector<double> vec1, vector<double> vec2){
    vector<double> output;
    for(int i=0;i<vec1.size();i++){
        output.push_back(vec1[i]-vec2[i]);
    }
    return output;
}

vector<double> Subtract1DVectorsDoubInt(vector<double> vec1, vector<int> vec2){
    vector<double> output;
    for(int i=0;i<vec1.size();i++){
        output.push_back(vec1[i]-(double)vec2[i]);
    }
    return output;
}

vector<int> BinaryPSF(vector<vector<double>> array_of_PSFs, vector<int> secondary_coordinates, vector<double> image_psf, double flux, vector<double> w)
{
    //cout << secondary_coordinates << endl;
	vector<double> scale_array(19);
    for(int i= 0.95; i > 0; i = i-0.05){
        scale_array.push_back(i);
    }
	int models_length = array_of_PSFs.size();
	int scale_length = scale_array.size();
	int chi_array_length = models_length * models_length * scale_length;
	
	vector<double> chi_array;
	chi_array.reserve(chi_array_length);
	
	int cx = (array_of_PSFs.size() - 1) * .5;
	int cy = (array_of_PSFs[0].size() - 1) * .5;
	cout << "line 66" << endl;
	for (int i = 0; i < models_length; i++) { 
		vector<double> model_PSF = array_of_PSFs[i];
        vector<double> cx_range;
        vector<double> cy_range;
        vector<double> PSF1;
        for (int j = cy-2; j < cy+3; j++){
            PSF1.push_back(j);
        }
        for (int j = cx-2; j < cx+3; j++){
            PSF1.push_back(j);
        }
        //PSF1.push_back(cy_range);
        //PSF1.push_back(cx_range);
        //A slice of model_PSF from model_PSF[cy-2:cy+3,cx-2:cx+3]
		//vector<double> PSF1 = (xt::view(model_PSF, xt::range(cy - 2, cy + 3 ), xt::range(cx - 2, cx + 3)));
		for (int j = 0; j < models_length; ++j) {
			vector<double> model_PSF = array_of_PSFs[j];
			int rangeY_1 = cy - secondary_coordinates[0];
			int rangeY_2 = cy - secondary_coordinates[0] + 5;
			int rangeX_1 = cx - secondary_coordinates[1];
			int rangeX_2 = cx - secondary_coordinates[1] + 5;
			
            vector<double> cx_range2;
            vector<double> cy_range2;
            vector<double> PSF2;
            for (int j = rangeY_1; j < rangeY_2; j++){
                PSF2.push_back(j);
            }
            for (int j = rangeX_1; j < rangeX_2; j++){
                PSF2.push_back(j);
            }
            //PSF2.push_back(cy_range2);
            //PSF2.push_back(cx_range2);
			//vector<double> PSF2 = (xt::view(model_PSF, xt::range(rangeY_1, rangeY_2), xt::range(rangeX_1, rangeX_2)));	
            
			for (int k = 0; k < scale_length; ++k) {          
				double relative_flux = scale_array[k];
                
                vector<double> a = PSF1;
                vector<double> b = PSF2;
                a = MultiplyByScalar1D(a, relative_flux);
                b = MultiplyByScalar1D(b, (1-relative_flux));
                vector<double> PSF_sum = Add1DVectors(a,b);
                vector<double> temp_array = MultiplyByScalar1D(PSF_sum,flux);
                
                vector<double> chi_PSF = Subtract1DVectors(image_psf,temp_array);
                chi_PSF = Multiply1DVectors(chi_PSF, w);
				//vector<double> chi_PSF = w * (image_psf - temp_array);
                double chi_value = SumSquared1D(chi_PSF);
				//double chi_value = sum(chi_PSF * chi_PSF)(0);
				chi_array.push_back(chi_value);	
			}
		}
	} 
    cout << "line 89" << endl;
    vector<vector<vector<double>>> chi = Reshape3D(chi_array, 100,100,19);
	//auto chi = xt::adapt(chi_array, {100, 100, 19});
    //This is a reshape function, reshaping chi's sorted PSFs by best to worst psf's
    vector<int> best_PSFs{0,0,0};
    for (int i=0; i<chi.size();i++){
        for (int j=0; j<chi[i].size();j++){
            for (int k=0; k<chi[i][j].size();k++){
                double new_val = chi[i][j][k];
                double best = chi[best_PSFs[i]][best_PSFs[j]][best_PSFs[k]];
                if (new_val < best) { 
                    best_PSFs[0] = i;
                    best_PSFs[1] = j;
                    best_PSFs[2] = k;
                }
            }
        }
    }
	//vector<int> best_PSFs = xt::adapt(xt::unravel_index(xt::argmin(chi)(0), chi.size(), xt::layout_type::row_major), {3});
	//cout << "line 91: " << best_PSFs << endl;
	return best_PSFs;
}

vector<double> BinaryRelativeFlux(vector<vector<double>> array_of_PSFs, vector<int> secondary_coordinates, vector<double> image_psf, double flux, vector<int> w, int best_primary, int best_secondary) 
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
    for(int i = 0.99; i > 0; i=i-0.01){
        scale_array.push_back(i);
    }
	//vector<double> scale_array = xt::arange<double>(.99, 0, -.01);
	
	int primary_length = primary_ref.size();
	int secondary_length = secondary_ref.size();
	int scale_length = scale_array.size();
	int chi_array_length = primary_length * secondary_length * scale_length;
	vector<double> chi_array;
	chi_array.reserve(chi_array_length);
	
	for (int i = 0; i < primary_length; ++i ) {
		int primary = primary_ref[i];
		vector<double> model_PSF = array_of_PSFs[primary];
		int cx = (model_PSF.size() - 1) * .5;
		int cy = (model_PSF.size() - 1) * .5;
        vector<double> cx_range;
        vector<double> cy_range;
        vector<double> PSF1;
        for (int j = cy-2; j < cy+3; j++){
            PSF1.push_back(j);
        }
        for (int j = cx-2; j < cx+3; j++){
            PSF1.push_back(j);
        }
        //PSF1.push_back(cy_range);
        //PSF1.push_back(cx_range);
		//vector<double> PSF1 = (xt::view(model_PSF, xt::range(cy - 2, cy + 3 ), xt::range(cx - 2, cx + 3)));	
		for (int j = 0; j < secondary_length; ++j) {
			int secondary = secondary_ref[j];
			vector<double> model_PSF = array_of_PSFs[secondary];
			int rangeY_1 = cy - secondary_coordinates[0];
			int rangeY_2 = cy - secondary_coordinates[0] + 5;
			int rangeX_1 = cx - secondary_coordinates[1];
			int rangeX_2 = cx - secondary_coordinates[1] + 5;
			
            vector<double> cx_range2;
            vector<double> cy_range2;
            vector<double> PSF2;
            for (int j = rangeY_1; j < rangeY_2; j++){
                PSF2.push_back(j);
            }
            for (int j = rangeX_1; j < rangeX_2; j++){
                PSF2.push_back(j);
            }
            //PSF2.push_back(cy_range2);
            //PSF2.push_back(cx_range2);
			//vector<double> PSF2 = (xt::view(model_PSF, xt::range(rangeY_1, rangeY_2), xt::range(rangeX_1, rangeX_2)));	
			for (int k = 0; k < scale_length; ++k) {
				double relative_flux = scale_array[k];
				vector<double> a = MultiplyByScalar1D(PSF1, relative_flux);
				vector<double> b = MultiplyByScalar1D(PSF2, (1 - relative_flux));
				vector<double> temp_array = MultiplyByScalar1D(Add1DVectors(a, b), flux);
				vector<double> chi_PSF = Multiply1DVectorsInt(w, Subtract1DVectors(image_psf , temp_array));
                //double chi_value;
                //for(int i = 0; i < chi_PSF.size(), i++){
                    
                //}
                
				double chi_value = SumSquared1D(chi_PSF);
				chi_array.push_back(chi_value);
			}
		}
	}
	
	//vector<double> chi = xt::adapt(chi_array, {9, 9, 99});
    vector<vector<vector<double>>> chi = Reshape3D(chi_array, 9, 9, 99);
    vector<int> best_PSFs = FindBestPSF(chi);
	//vector<int> best_PSFs = xt::adapt(xt::unravel_index(xt::argmin(chi)(0), chi.size(), xt::layout_type::row_major), {3});
	double best_primary_ref = primary_ref[best_PSFs[0]];
	double best_secondary_ref = secondary_ref[best_PSFs[1]];
	double best_ref = scale_array[best_PSFs[2]];
	vector<double> final_best_PSF = {best_primary_ref, best_secondary_ref, best_ref};
	
	return final_best_PSF;
}

double BinaryFit(vector<double> model_primary, vector<double> model_secondary, vector<int> secondary_coordinates, vector<double> image_psf, double base, vector<int> w, double relative_flux) 
{
	int cx = (model_primary.size() - 1) * .5;
	int cy = (model_primary.size() - 1) * .5;
	int rangeY_1 = cy - secondary_coordinates[0];
	int rangeY_2 = cy - secondary_coordinates[0] + 5;
	int rangeX_1 = cx - secondary_coordinates[1];
	int rangeX_2 = cx - secondary_coordinates[1] + 5;
	
    vector<double> PSF1;
    for (int j = cy-2; j < cy+3; j++){
        PSF1.push_back(j);
    }
    for (int j = cx-2; j < cx+3; j++){
        PSF1.push_back(j);
    }
    //PSF1.push_back(cy_range);
    //PSF1.push_back(cx_range);
    vector<double> PSF2;
    for (int j = rangeY_1; j < rangeY_2; j++){
        PSF2.push_back(j);
    }
    for (int j = rangeX_1; j < rangeX_2; j++){
        PSF2.push_back(j);
    }
    //PSF2.push_back(cy_range);
    //PSF2.push_back(cx_range);
	//vector<double> PSF1 = (xt::view(model_primary, xt::range(cy - 2, cy + 3 ), xt::range(cx - 2, cx + 3)));		
	//vector<double> PSF2 = (xt::view(model_secondary, xt::range(rangeY_1, rangeY_2), xt::range(rangeX_1, rangeX_2)));
	vector<double> flux_array;
	vector<double> chi_array;
	
	//change "range" to something like "Scale"
    
	vector<double> scale_array(1001); // = xt::arange<double>(0, 10.01, 0.01);
    int j=0;
    for (int i=0;i<10.01;i=i+0.01){
        scale_array[j] = i;
        j++;
    }
	int scale_length = scale_array.size();
	flux_array.reserve(scale_length);
	chi_array.reserve(scale_length);
	vector<double> PSF_sum = Add1DVectors(MultiplyByScalar1D(PSF1, relative_flux), MultiplyByScalar1D(PSF2, (1 - relative_flux)));
	double scale;
	for (int i = 0; i < scale_length; ++i) {
		scale = scale_array[i];
		vector<double> scale_PSF = MultiplyByScalar1D(MultiplyByScalar1D(PSF_sum, scale), base);
		flux_array.push_back(scale * base);
		vector<double> chi_PSF = Multiply1DVectorsInt(w, Subtract1DVectors(image_psf, scale_PSF));
		double chi_value = SumSquared1D(chi_PSF);
		chi_array.push_back(chi_value);
	}
		
	auto chi_squareds = xt::adapt(chi_array, {scale_length});
	auto fluxes = xt::adapt(flux_array, {scale_length});
    
	double min_chi = chi_squareds(0);
    for(int i=1;i<chi_squareds.size();i++){
        if (chi_squareds(i) < min_chi) { min_chi = chi_squareds(i); }
    }
	int min_chi_index = 0;
	int chi_length = chi_squareds.size();
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

double SecondaryBinaryFit(vector<double> primary, vector<double> secondary, vector<int> secondary_coordinates, vector<int> w2, vector<double> image_psf, double flux, double relative_flux)
{
	int cx = (primary.size() - 1) * .5;
	int cy = (primary.size() - 1) * .5;
	int rangeY_1 = cy - secondary_coordinates[0];
	int rangeY_2 = cy - secondary_coordinates[0] + 5;
	int rangeX_1 = cx - secondary_coordinates[1];
	int rangeX_2 = cx - secondary_coordinates[1] + 5;
	vector<double> PSF1;
    vector<double> PSF2;
    for (int j = cy-2; j < cy+3; j++){
        PSF1.push_back(j);
    }
    for (int j = cx-2; j < cx+3; j++){
        PSF1.push_back(j);
    }
    //PSF1.push_back(cy_range);
    //PSF1.push_back(cx_range);
    
    for (int j = rangeY_1; j < rangeY_2; j++){
        PSF2.push_back(j);
    }
    for (int j = rangeX_1; j < rangeX_2; j++){
        PSF2.push_back(j);
    }
    //PSF2.push_back(cy_range);
    //PSF2.push_back(cx_range);
	//vector<double> PSF1 = (xt::view(primary, xt::range(cy - 2, cy + 3), xt::range(cx - 2, cx + 3)));
	//vector<double> PSF2 = (xt::view(secondary, xt::range(rangeY_1, rangeY_2), xt::range(rangeX_1, rangeX_2)));
	vector<double> scale_array(100001);// = xt::arange<double>(0, 10.0001, .0001);
    int j = 0;
    for(int i=0;i<10.0001;i=i+0.0001){
        scale_array[j] = i;
    }
    
    
	int scale_length = scale_array.size();
	vector<double> flux_array;
	vector<double> chi_array;
	flux_array.reserve(scale_length);
	chi_array.reserve(scale_length);
    
	
	vector<double> relative_PSF2 = Multiply1DVectorsInt(w2, Subtract1DVectors(image_psf, MultiplyByScalar1D(PSF1, flux*relative_flux)));
                                                     
                                                     
	for (int i = 0; i < scale_length; ++i)
	{
		double scale = scale_array[i];
		vector<double> scale_PSF2 =  MultiplyByScalar1D(MultiplyByScalar1D(Multiply1DVectorsInt(w2,PSF2),scale),flux);
		flux_array.push_back(flux* scale);
		vector<double> chi_PSF = Multiply1DVectorsInt(w2,Subtract1DVectors(relative_PSF2,scale_PSF2));
		double chi_value = SumSquared1D(chi_PSF);
		chi_array.push_back(chi_value);
	}

	auto chi_squareds = xt::adapt(chi_array, {scale_length});
	auto fluxes = xt::adapt(flux_array, {scale_length});
    double min_chi = chi_squareds(0);
    for(int i=1;i<chi_squareds.size();i++){
        if (chi_squareds(i) < min_chi) { min_chi = chi_squareds(i); }
    }
	//double min_chi = xt::amin(chi_squareds)(0);
	int min_chi_index = 0;
	int chi_length = chi_squareds.size();
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
	int cx = (PSF_1.size() - 1) * .5;
	int cy = (PSF_2.size() - 1) * .5;
	
	int P_rangeY_1 = cy - coordinates_1[0];
	int P_rangeY_2 = cy - coordinates_1[0] + 5;
	int P_rangeX_1 = cx - coordinates_1[1];
	int P_rangeX_2 = cx - coordinates_1[1] + 5;
	
	int S_rangeY_1 = cy - coordinates_2[0];
	int S_rangeY_2 = cy - coordinates_2[0] + 5;
	int S_rangeX_1 = cx - coordinates_2[1];
	int S_rangeX_2 = cx - coordinates_2[1] + 5;

    vector<int> PSF1;
    vector<double> cy_range;
    vector<double> cx_range;
    for (int j = P_rangeY_1; j < P_rangeY_2; j++){
        PSF1.push_back(j);
    }
    for (int j = P_rangeX_1; j < P_rangeX_2; j++){
        PSF1.push_back(j);
    }
    //PSF1.push_back(cy_range);
    //PSF1.push_back(cx_range);
    vector<int> PSF2;
    for (int j = S_rangeY_1; j < S_rangeY_2; j++){
        PSF2.push_back(j);
    }
    for (int j = S_rangeX_1; j < S_rangeX_2; j++){
        PSF2.push_back(j);
    }
    //PSF2.push_back(cy_range);
    //PSF2.push_back(cx_range);
	//vector<double> PSF1 = (xt::view(PSF_1, xt::range(P_rangeY_1, P_rangeY_2), xt::range(P_rangeX_1, P_rangeX_2)));
	//vector<double> PSF2 = (xt::view(PSF_2, xt::range(S_rangeY_1, S_rangeY_2), xt::range(S_rangeX_1, S_rangeX_2)));
	
	vector<double> scale_array(100001);// = xt::arange<double>(0, 10.0001, .0001);
    int j = 0;
    for(int i =0;i<10.001;i++){
        scale_array[j]= i;
    }
	int scale_length = scale_array.size();
	vector<double> flux_array;
	vector<double> chi_array;
	flux_array.reserve(scale_length);
	chi_array.reserve(scale_length);
	
	vector<double> rPSF = Multiply1DVectorsInt(w, Subtract1DVectorsDoubInt(image_psf, MultiplyByScalar1DInt(PSF2, flux)));
	for (int i = 0; i < scale_length; ++i) {
		double scale = scale_array[i];
		vector<double> scale_PSF = MultiplyByScalar1D(MultiplyByScalar1D(Multiply1DVectorsIntInt(PSF1,w),scale),base);
		flux_array.push_back(scale * base);
		vector<double> chi_PSF = Multiply1DVectorsInt(w, Subtract1DVectors(rPSF,scale_PSF));
		double chi_value = SumSquared1D(chi_PSF);
		chi_array.push_back(chi_value);
	}
	
	auto chi_squareds = xt::adapt(chi_array, {scale_length});
	auto fluxes = xt::adapt(flux_array, {scale_length});
	//double min_chi = xt::amin(chi_squareds)(0);
    double min_chi = chi_squareds[0];
    for(int i=1;i<chi_squareds.size();i++){
        if (chi_squareds[i] < min_chi) { min_chi = chi_squareds[i]; }
    }
	int min_chi_index = 0;
	int chi_length = chi_squareds.size();
	for (int i = 0; i < chi_length ; ++i) {
		if (chi_squareds[i] == min_chi) {
			min_chi_index = i;
			break;
		}
	}	
	
	double best_flux = fluxes(min_chi_index);
	
	return best_flux;
}

vector<vector<int>> Read_FileInt(string filename){
    vector<vector<int>> vec; 
    short loop=0; 
    //vector<double> line; 
    ifstream myfile(filename); 
    int row, col;
    if (myfile.is_open()) 
    {
        myfile >> row >> col;
        for(int i = 0; i < row; i++){
            vector<int> temp;
            for(int j = 0; j < col; j++ ){
                int temp_val;
                myfile >> temp_val;
                temp.push_back(temp_val);
                //cout << vec[i][j];
            }
            vec.push_back(temp);
        }
        myfile.close();
    }
    return vec;
}
         

vector<vector<double>> Read_FileDouble(string filename){
    vector<vector<double>> vec; 
    short loop=0; 
    //vector<double> line; 
    ifstream myfile(filename); 
    int row, col;
    if (myfile.is_open()) 
    {
        myfile >> row >> col;
        for(int i = 0; i < row; i++){
            vector<double> temp;
            for(int j = 0; j < col; j++ ){
                double temp_val;
                myfile >> temp_val;
                temp.push_back(temp_val);
                cout << vec[i][j];
            }
            vec.push_back(temp);
        }
        myfile.close();
    }
    return vec;
}

vector<double> Read_FileDouble1D(string filename){
    vector<double> vec; 
    short loop=0; 
    //vector<double> line; 
    ifstream myfile(filename); 
    int row;
    if (myfile.is_open()) 
    {
        myfile >> row;
        for(int i = 0; i < row; i++){
            double temp_val;
            myfile >> temp_val;
            vec.push_back(temp_val);
            cout << vec[i];
        }
        //vec.push_back(temp);    
        myfile.close();
    }
    return vec;
}


//MAIN 
int main() {    
    //read in the image and model data, and the parameters we set earlier
    string ipsf = "tempdata/image_psf.txt";
    vector<double> image_psf = Read_FileDouble1D(ipsf);   
    cout << "Secondary line 820 is working!" << endl;
    int array_length = 100;
    vector<vector<double>> array_of_PSFs(100);
    //double array_of_PSFs = xt::load_npy<double>("tempdata/array_of_PSFs.txt");
    	for (int i = 1; i <= array_length; ++i) {
    		string filename = "tempdata/PSFmodel_" + to_string(i) + ".txt";
    		array_of_PSFs[i - 1] = Read_FileDouble1D(filename);
    	}
    	
    string seccoords = "tempdata/sec_coords.txt";
    string primcoords = "tempdata/prim_coords.txt";
	vector<vector<int>> coordinates = Read_FileInt(seccoords);
	vector<vector<double>> center_load = Read_FileDouble(primcoords);
    cout << center_load[0][0] << endl;
	vector<int> center = {(int)center_load[0][0], (int)center_load[1][1]};
    	
    	vector<int> primary_coordinates = {2,2};
    	vector<int> secondary_coordinates = vector<int>(2);
    	
    	//arrays to store the output values
    	vector<double> residual_error_array(coordinates.size());
	vector<int> best_primary_array(coordinates.size());
	vector<double> flux_primary_array(coordinates.size());
	vector<int> best_secondary_array(coordinates.size());
	vector<double> flux_secondary_array(coordinates.size());
	vector<double> flux_sum_array(coordinates.size());
    	
    
    	//do the binary fit
    	cout << "Beginning Second Binary Fit" << endl;
        //cout << "coordinates" << coordinates << "center" << center << endl;
    
    //The coordinates provided are the 9 pixels surrounding the centroid
    //For the second binary fit, we want to provide the 9 subpixels that surround the best fit primary and the secondary, in this case, 141.73,150.45 and 143.23,150.57
    	for (int i = 0; i < coordinates.size(); ++i) {
            cout << "loop " << i << endl;
    		vector<int> secondary = coordinates[i];
    		secondary_coordinates[0] = secondary[0] - center[0] + 2;
    		secondary_coordinates[1] = secondary[1] - center[1] + 2;
    		//Create a new Binary Base from the new coordinates given
    		vector<vector<int>> w_array = BinaryBase(primary_coordinates, secondary_coordinates);
            
            vector<int> w = w_array[0];       
            vector<int> w1 = w_array[1];       
            vector<int> w2 = w_array[2];       
            
		//vector<int> w1 = vector<int> dest(w_array[1].begin(), w_array[1].end());
		//vector<int> w2 = vector<int> dest(w_array[2].begin(), w_array[2].end());
            vector<double> sqrt_weights;
            vector<int> inted_sqrtw;
            for(int i=0;i<image_psf.size();i++){
                double val = image_psf[i];
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
    		double flux = base * 1;
    		cout << "line 395" << endl;
    		while (abs(comp - flux) > .000001 && iteration < 10) {
    			comp = flux;
                cout << "line 398" << endl;
    			vector<int> best_PSFs =  BinaryPSF(array_of_PSFs, secondary_coordinates, image_psf, flux, sqrt_weights);
                cout << "line 401" << endl;
    			best_primary = best_PSFs[0];
    			best_secondary = best_PSFs[1];
                cout << "line 402" << endl;
                //vector<int> inted_sqrtw = vector<int> dest(sqrt_weights.begin(), sqrt_weights.end());
    			vector<double> best_PSFs_ref = BinaryRelativeFlux(array_of_PSFs, secondary_coordinates, image_psf, flux, inted_sqrtw, best_primary, best_secondary);
    			best_primary = best_PSFs_ref[0];
    			best_secondary = best_PSFs_ref[1];
    			best_RF = best_PSFs_ref[2];
    			cout << "line 407" << endl;
			flux = BinaryFit(array_of_PSFs[best_primary], array_of_PSFs[best_secondary], secondary_coordinates, image_psf, base, w, best_RF);
                cout << "line 409" << endl;
    			++iteration;
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
                if (iteration >= 10) {break;}
			//CREATE rfpsf FUNCTION
                flux_primary = RelativeFluxes(primary_PSF, primary_coordinates, secondary_PSF, secondary_coordinates, w1, image_psf, base, flux_secondary);
                flux_secondary = RelativeFluxes(secondary_PSF, secondary_coordinates, primary_PSF, primary_coordinates, w2, image_psf, base, flux_primary);
            
                ++iteration;
            }            
            
		
		int cx = (primary_PSF.size() - 1) * .5;
		int cy = (primary_PSF.size() - 1) * .5;
		
		int P_rangeY_1 = cy - primary_coordinates[0];
		int P_rangeY_2 = cy - primary_coordinates[0] + 5;
		int P_rangeX_1 = cx - primary_coordinates[1];
		int P_rangeX_2 = cx - primary_coordinates[1] + 5;
		
		int S_rangeY_1 = cy - secondary_coordinates[0];
		int S_rangeY_2 = cy - secondary_coordinates[0] + 5;
		int S_rangeX_1 = cx - secondary_coordinates[1];
		int S_rangeX_2 = cx - secondary_coordinates[1] + 5;
            
       
        vector<double> PSF1;
        //PSF1.reserve(P_rangeY_2-P_rangeY_1);
        vector<double> PSF2;
        //PSF2.reserve(S_rangeY_2-S_rangeY_1);
        for (int i = P_rangeY_1; i < P_rangeY_2; ++i) {
            PSF1.push_back(primary_PSF[i]);
            //PSF1.emplace_back(primary_PSF[i].begin() + P_rangeX_1, primary_PSF[i].begin() + P_rangeX_2);
        }
        for (int i = S_rangeY_1; i < S_rangeY_2; ++i) {
            PSF2.push_back(primary_PSF[i]);
            //PSF2.emplace_back(primary_PSF[i].begin() + S_rangeX_1, primary_PSF[i].begin() + S_rangeX_2);
        }
            
		//vector<double> PSF1 = slicing(primary_PSF, xt::range(P_rangeY_1, P_rangeY_2), xt::range(P_rangeX_1, P_rangeX_2)));
		//vector<double> PSF2 = (xt::view(secondary_PSF, xt::range(S_rangeY_1, S_rangeY_2), xt::range(S_rangeX_1, S_rangeX_2)));

		vector<double> PSF_sum = Add1DVectors(MultiplyByScalar1D(PSF1,flux_primary), MultiplyByScalar1D(PSF2,flux_secondary));
            
		vector<double> pre_residual = Subtract1DVectors(image_psf,PSF_sum);
        for(int i=0;i< pre_residual.size();i++){
            if(pre_residual[i] < 0) { pre_residual[i] = pre_residual[i]*-1; }
        }
        vector<double> temp_vec = Multiply1DVectorsInt(w,pre_residual);
		double residual = SumVector(temp_vec);
		double residual_error = residual / SumVectorInt(w);
		
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
    	vector<vector<int>> center_array = {{center[0], center[1]}, {center[0], center[1]}, {center[0], center[1]}, {center[0], center[1]}, {center[0], center[1]}, {center[0], center[1]}, {center[0], center[1]}, {center[0], center[1]}, {center[0], center[1]}};
    	vector<vector<int>> secondary_array = coordinates;
    	
    	//send the "outputs" as seperate files since they are different data types and we can read them all in to python for the calculation of angles and separation - it only takes a second :)
    ofstream outFile1("tempdata/residual_error_array.txt");
    for (const auto &e : residual_error_array) outFile1 << e << "\n";
    ofstream outFile2("tempdata/center_array.txt");
    for (const vector<int> &e : center_array){ 
        for (const int &f : e){
            outFile2 << f << " ";
        }
        outFile2 << "/n";
    }
    ofstream outFile3("tempdata/secondary_array.txt");
    for (const vector<int> &e : secondary_array){ 
        for (const int &f : e){
            outFile3 << f << " ";
        }
        outFile3 << "/n";
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
