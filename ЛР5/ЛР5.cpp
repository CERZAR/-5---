#include "pch.h"
#include <numeric>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <exception>
#include <map>

#define line 	std::cout << std::endl;

void string_delimetr(int &&max, char ch)
{
	for (size_t i = 1; i < max; i++)
		std::cout << ch;
	line;
}

struct alternative
{
	std::string name;
	int appearance = NULL;
	int cost = NULL;
	int home = NULL;
	int character = NULL;
};

struct criterion
{
	std::string name;
	bool is_main = false;
};

void init1(std::vector<alternative> &alternatives)
{
	alternative a;
	a.name = "Tatiana";
	a.appearance = 6;
	a.cost = 2;
	a.home = 3;
	a.character = 4;
	alternatives.push_back(a);
	a.name = "Larisa";
	a.appearance = 8;
	a.cost = 7;
	a.home = 2;
	a.character = 8;
	alternatives.push_back(a);
	a.name = "Natalia";
	a.appearance = 5;
	a.cost = 3;
	a.home = 9;
	a.character = 6;
	alternatives.push_back(a);
	a.name = "Olga";
	a.appearance = 7;
	a.cost = 2;
	a.home = 8;
	a.character = 2;
	alternatives.push_back(a);
}

void init2(std::vector<criterion> &criteria)
{
	criterion a;
	a.name = "Appearance";
	criteria.push_back(a);
	a.name = "Expenses";
	criteria.push_back(a);
	a.name = "Housekeep";
	criteria.push_back(a);
	a.name = "Character";
	criteria.push_back(a);
}

void print_vec_weight_criteria(std::vector<criterion> criteria, std::vector <double> vec_weight_criteria) // table 1
{
	if (criteria.size() != vec_weight_criteria.size())
		throw std::exception("Vectors size error");
	std::cout << "Vector of weights of criteria:";
	line;
	string_delimetr(42, '_');
	std::cout << "Criteria	Weight	Normalized_vector";
	line;
	double sum = accumulate(vec_weight_criteria.begin(),  vec_weight_criteria.end(), 0);
	for (size_t i = 0; i < criteria.size(); ++i)
		std::cout << criteria[i].name << ":	  " << vec_weight_criteria[i] << "		" << vec_weight_criteria[i] / sum << std::endl;
	string_delimetr(42, '_');
	line;
}

int method_1(std::vector<alternative> alternatives, std::vector<criterion> criteria, std::vector <double> vec_weight_criteria)
{
	std::cout << "	#1) Method for replacing criteria by limitations:\n\n\n";
	const int main_criterion = 3;
	criteria[main_criterion].is_main = true;
	std::cout << "The main criterion was chosen - " << criteria[main_criterion].name;
	line;
	std::cout << "Matrix of alternatives assessment:\n";
	string_delimetr(71, '_');
	for (size_t i = 0; i < alternatives.size(); ++i)
	{
		std::cout << "	    " << criteria[i].name;
		if (criteria[i].is_main == true)
			std::cout << "!";
	}
	line;
	for (size_t i = 0; i < alternatives.size(); ++i)
		std::cout << alternatives[i].name << "		" << alternatives[i].appearance << "		" << alternatives[i].cost 
		<< "		" << alternatives[i].home << "		" << alternatives[i].character << std::endl;
	string_delimetr(71, '_');
	std::vector<std::vector<int>> criteria_vector;
	std::vector<double> min_conditions;
	std::vector<int> a;
		for (size_t j = 0; j < alternatives.size(); ++j)
			a.push_back(alternatives[j].appearance);
		criteria_vector.push_back(a);
		a.clear();
		min_conditions.push_back(vec_weight_criteria[0] * 0.1 * *std::max_element(criteria_vector[0].begin(), criteria_vector[0].end()));
		for (size_t j = 0; j < alternatives.size(); ++j)
			a.push_back(alternatives[j].cost);
		criteria_vector.push_back(a);
		a.clear();
		min_conditions.push_back(vec_weight_criteria[1] * 0.1 * *std::max_element(criteria_vector[1].begin(), criteria_vector[1].end()));
		for (size_t j = 0; j < alternatives.size(); ++j)
			a.push_back(alternatives[j].home);
		criteria_vector.push_back(a);
		a.clear();
		min_conditions.push_back(vec_weight_criteria[2] * 0.1 * *std::max_element(criteria_vector[2].begin(), criteria_vector[2].end()));
		for (size_t j = 0; j < alternatives.size(); ++j)
			a.push_back(alternatives[j].character);
		criteria_vector.push_back(a);
		min_conditions.push_back(vec_weight_criteria[3] * 0.1 * *std::max_element(criteria_vector[3].begin(), criteria_vector[3].end()));
	line;
	std::cout << "Minimum acceptable values for other criteria:\n";
	string_delimetr(46, '_');
	for (size_t i = 0; i < alternatives.size(); ++i)
		if (criteria[i].is_main != true)
			std::cout << criteria[i].name << " no less than: " << vec_weight_criteria[i] * 0.1 << " * " 
			<< *std::max_element(criteria_vector[i].begin(), criteria_vector[i].end()) << " = " << min_conditions[i] << std::endl;
	string_delimetr(46, '_');
	std::vector<int>current_candidates({ 0,1,2,3 });
	for (size_t i = 0; i < alternatives.size(); ++i)
	{
		if ((criteria[0].is_main != true) && (alternatives[i].appearance < min_conditions[0]))
			current_candidates.erase(std::find(current_candidates.begin(), current_candidates.end(), i));
		else if ((criteria[1].is_main != true) && (alternatives[i].cost < min_conditions[1]))
			current_candidates.erase(std::find(current_candidates.begin(), current_candidates.end(), i));
		else if ((criteria[2].is_main != true) && (alternatives[i].home < min_conditions[2]))
			current_candidates.erase(std::find(current_candidates.begin(), current_candidates.end(), i));
		else if ((criteria[3].is_main != true) && (alternatives[i].character < min_conditions[3]))
			current_candidates.erase(std::find(current_candidates.begin(), current_candidates.end(), i));
	}
	if (current_candidates.size() > 1)
	{
		std::vector<int> sum;
		std::vector<int> result_candidates;
		for (size_t i = 0; i < current_candidates.size(); ++i)
		{
			int s = 0;
			if (criteria[0].is_main != true)
				s += alternatives[current_candidates[i]].appearance;
			if (criteria[1].is_main != true)
				s += alternatives[current_candidates[i]].cost;
			if (criteria[2].is_main != true)
				s += alternatives[current_candidates[i]].home;
			if (criteria[3].is_main != true)
				s += alternatives[current_candidates[i]].character;
			sum.push_back(s);
		}
		int max = *std::max_element(sum.begin(), sum.end());
		int counter = 0;
		for (size_t i = 0; i < sum.size(); ++i)
		{
			if (sum[i] == max)
			{
				counter++;
				result_candidates.push_back(current_candidates[i]);
			}
		}
		if (counter > 1)
		{
			sum.clear();
			for (size_t i = 0; i < result_candidates.size(); ++i)
			{
				if (criteria[0].is_main == true)
					sum.push_back(alternatives[result_candidates[i]].appearance);
				if (criteria[1].is_main == true)
					sum.push_back(alternatives[result_candidates[i]].cost);
				if (criteria[2].is_main == true)
					sum.push_back(alternatives[result_candidates[i]].home);
				if (criteria[3].is_main == true)
					sum.push_back(alternatives[result_candidates[i]].character);
			}
			max = *std::max_element(sum.begin(), sum.end());
			for (size_t i = 0; i < sum.size(); ++i)
			{
				if (sum[i] == max)
				{
					std::cout << "Found several suitable alternatives: ";
					for (size_t i = 0; i < result_candidates.size(); ++i)
					{
						std::cout << alternatives[result_candidates[i]].name;
						if (i != result_candidates.size() - 1)
							std::cout << ", ";
						else
							std::cout << ". ";
					}
					std::cout << " Checking by main criterion...\n";
					std::cout << "\n|> The most suitable alternative - " << alternatives[result_candidates[i]].name;
				}
			}
		}
		else
			std::cout << "\n|> The only suitable alternative - " << alternatives[result_candidates[0]].name;
	}
	else if (current_candidates.size() == 0)
	{
		int max = 0;
		std::cout << "No suitable alternatives found. Checking by main criterion...\n";
		if (criteria[0].is_main == true)
		{
			max = *std::max_element(criteria_vector[0].begin(), criteria_vector[0].end());
			for (size_t i = 0; i < alternatives.size(); ++i)
				if (alternatives[i].appearance == max)
					std::cout << "\n|> The only suitable alternative - " << alternatives[i].name;
		}
		if (criteria[1].is_main == true)
		{
			max = *std::max_element(criteria_vector[1].begin(), criteria_vector[1].end());
			for (size_t i = 0; i < alternatives.size(); ++i)
				if (alternatives[i].cost == max)
					std::cout << "\n|> The only suitable alternative - " << alternatives[i].name;
		}
		if (criteria[2].is_main == true)
		{
			max = *std::max_element(criteria_vector[2].begin(), criteria_vector[2].end());
			for (size_t i = 0; i < alternatives.size(); ++i)
				if (alternatives[i].home == max)
					std::cout << "\n|> The only suitable alternative - " << alternatives[i].name;
		}
		if (criteria[3].is_main == true)
		{
			max = *std::max_element(criteria_vector[3].begin(), criteria_vector[3].end());
			for (size_t i = 0; i < alternatives.size(); ++i)
				if (alternatives[i].character == max)
					std::cout << "\n|> The only suitable alternative - " << alternatives[i].name;
		}
	}
	else
		std::cout << "\n|> The only suitable alternative - " << alternatives[current_candidates[0]].name;
	criteria[main_criterion].is_main = false;
	return current_candidates[0];
}

double find_distance(double x, double y)
{
	return (sqrt(pow(10 - x, 2) + pow(10 - y, 2)));
}

int method_2(std::vector<alternative> alternatives, std::vector<criterion> criteria)
{
	std::cout << "\n\n\n	#2) The formation and narrowing of the Pareto set:\n\n\n";
	const int main_criterion_1 = 0;
	const int main_criterion_2 = 3;
	criteria[main_criterion_1].is_main = true;
	criteria[main_criterion_2].is_main = true;
	std::cout << "The main criteria was chosen - " << criteria[main_criterion_1].name << ", " << criteria[main_criterion_2].name;
	line;
	std::cout << "Utopia point - (10, 10)";
	line;
	std::cout << "Distance type - Euclidian";
	line;
	std::vector<double> distances;
	for (size_t i = 0; i < alternatives.size(); ++i)
		distances.push_back(find_distance(alternatives[i].appearance, alternatives[i].character));
	string_delimetr(53, '_');
	std::cout << "	    Appearence	    Character	    Distance" << std::endl;
	for (size_t i = 0; i < alternatives.size(); ++i)
	{
		std::cout << alternatives[i].name << "		" << alternatives[i].appearance << "		" << alternatives[i].character << "		" << std::fixed << std::setprecision(2) <<  distances[i] << std::endl;
	}
	string_delimetr(53, '_');
	int k = 0;
	double min = 11;
	for (size_t i = 0; i < distances.size()-1; ++i)
		if (distances[i] < min)
		{
			min = distances[i];
			k = i;
		}
	std::cout << "\n|> Alternative with minimal distance - " << alternatives[k].name << " (Distance: " << min << ")";
	return k;
}

std::vector<std::vector<int>> get_criteria_vector(std::vector<alternative> alternatives)
{
	std::vector<std::vector<int>> criteria_vector;
	std::vector<int> a;
	for (size_t j = 0; j < alternatives.size(); ++j)
		a.push_back(alternatives[j].appearance);
	criteria_vector.push_back(a);
	a.clear();
	for (size_t j = 0; j < alternatives.size(); ++j)
		a.push_back(alternatives[j].cost);
	criteria_vector.push_back(a);
	a.clear();
	for (size_t j = 0; j < alternatives.size(); ++j)
		a.push_back(alternatives[j].home);
	criteria_vector.push_back(a);
	a.clear();
	for (size_t j = 0; j < alternatives.size(); ++j)
		a.push_back(alternatives[j].character);
	criteria_vector.push_back(a);
	return criteria_vector;
}

int method_3(std::vector<alternative> alternatives, std::vector<criterion> criteria)
{
	std::cout << "\n\n\n	#3) Weighting and combining criteria:\n\n\n";
	std::vector<std::vector<int>> criteria_vector = get_criteria_vector(alternatives);
	std::vector<double> sum;
	for (size_t i = 0; i < criteria_vector.size(); ++i)
		sum.push_back(accumulate(criteria_vector[i].begin(), criteria_vector[i].end(), 0));
	std::cout << "Normalized matrix of alternatives assessment:  <-- matrix A" << std::endl;
	string_delimetr(71, '_');
	for (size_t i = 0; i < alternatives.size(); ++i)
		std::cout << "	     " << criteria[i].name;
	line;
	double normalized_alternatives_criteria_matrix[4][4];
	for (size_t i = 0; i < alternatives.size(); ++i)
	{
		normalized_alternatives_criteria_matrix[i][0] = alternatives[i].appearance / sum[0];
		normalized_alternatives_criteria_matrix[i][1] = alternatives[i].cost / sum[1];
		normalized_alternatives_criteria_matrix[i][2] = alternatives[i].home / sum[2];
		normalized_alternatives_criteria_matrix[i][3] = alternatives[i].character / sum[3];
	}
	for (size_t i = 0; i < alternatives.size(); ++i)
		std::cout << std::fixed << std::setprecision(4) << alternatives[i].name
		<< "	       " << normalized_alternatives_criteria_matrix[i][0]
		<< "	      " << normalized_alternatives_criteria_matrix[i][1]
		<< "	      " << normalized_alternatives_criteria_matrix[i][2]
		<< "	      " << normalized_alternatives_criteria_matrix[i][3]
		<< std::endl;
	double matrix[4][4];
	matrix[0][0] = 0.5;
	matrix[0][1] = 1;
	matrix[0][2] = 1;
	matrix[0][3] = 0.5;

	matrix[1][0] = 0;
	matrix[1][1] = 0.5;
	matrix[1][2] = 1;
	matrix[1][3] = 0;

	matrix[2][0] = 0;
	matrix[2][1] = 0;
	matrix[2][2] = 0.5;
	matrix[2][3] = 0;

	matrix[3][0] = 0;
	matrix[3][1] = 1;
	matrix[3][2] = 1;
	matrix[3][3] = 0.5;
	double vector[4];
	for (int i = 0; i < 4; ++i)
	{
		vector[i] = 0;
		for (int j = 0; j < 4; ++j)
			vector[i] += matrix[i][j];
		vector[i] -= 0.5;
	}
	double sumary = 0;
	for (int i = 0; i < 4; ++i)
		sumary += vector[i];
	string_delimetr(71, '_');
	std::cout << "Assessment of criteria: " << std::endl;
	string_delimetr(71, '_');
	for (size_t i = 0; i < alternatives.size(); ++i)
		std::cout << "	     " << criteria[i].name;
	line;
	std::cout << std::fixed << std::setprecision(1);
	for (size_t i = 0; i < alternatives.size(); ++i)
	{
		std::cout << criteria[i].name;
		for (size_t j = 0; j < 4; ++j)
		{
			std::cout << "	" << matrix[j][i] << "	";
		}
		line;
	}
	for (int i = 0; i < 4; ++i)
		vector[i] = vector[i] / sumary;
	string_delimetr(71, '_');
	std::cout << std::fixed << std::setprecision(2);
	std::cout << "Normalized criteria weight vector: ( ";
	for (int j = 0; j < 4; ++j)
		if (j != 3) std::cout << vector[j] << ", ";
		else std::cout << vector[j];
	std::cout << " )  <-- matrix B\n";
	double c[4];
	for (int i = 0; i < 4; i++)
	{
		c[i] = 0;
		for (int k = 0; k < 4; k++)
			c[i] += normalized_alternatives_criteria_matrix[i][k] * vector[k];
	}
	std::cout << "Multiply matrix A with matrix B...\n";

	std::cout << "\nAssesments for alternatives:\n";
	string_delimetr(71, '_');
	for (size_t i = 0; i < alternatives.size(); ++i)
		std::cout << "	     " << alternatives[i].name;
	line;
	std::cout << "Assesment";
	for (size_t i = 0; i < 4; ++i)
		std::cout << "     " << c[i] << "	 ";
	line;
	int k = 0;
	double max = 0;
	for (size_t i = 0; i < 4; ++i)
		if (c[i] > max)
		{
			max = c[i];
			k = i;
		}
	string_delimetr(71, '_');
	std::cout << "\n|> Alternative with the best assesment - " << alternatives[k].name << " (Assesment: " << max << ")";
	return k;
}

double get_mark(int delta)
{
	switch (delta)
	{
	case 0: return 1;
	case 1: return 3;
	case 2: return 4;
	case 3: return 5;
	case 4: return 6;
	case 5: return 7;
	case 6: return 8;
	case 7: return 9;
	case -1: return pow(3, -1);
	case -2: return pow(4, -1);
	case -3: return pow(5, -1);
	case -4: return pow(6, -1);
	case -5: return pow(7, -1);
	case -6: return pow(8, -1);
	case -7: return pow(9, -1);
	default: return 1;
	}
}

std::vector<double> normalise_matrix_by_string(std::vector<std::vector<double>> matrix)
{
	std::vector<double> vector;
	double sum = 0;
	for (size_t i = 0; i < 4; ++i)
	{
		for (size_t j = 0; j < 4; ++j)
			sum += matrix[i][j];
		vector.push_back(sum);
		sum = 0;
	}
	sum = 0;
	for (size_t i = 0; i < 4; ++i)
		sum += vector[i];
	for (size_t i = 0; i < 4; ++i)
		vector[i] = vector[i] / sum;
	return vector;
}

std::vector<double> multiply_matrix(std::vector<std::vector<double>> a, std::vector<double> b)
{
	std::vector<double> c;
	double current = 0;
	for (int i = 0; i < 4; i++)
	{
		for (int k = 0; k < 4; k++)
			current += a[i][k] * b[k];
		c.push_back(current);
		current = 0;
	}
	return c;
}

int method_4(std::vector<alternative> alternatives, std::vector<criterion> criteria, std::vector<double> vec_weight_criteria)
{
	std::cout << "\n\n\n	#4) Hierarchy Analysis Method:\n\n\n";
	std::vector<std::vector<int>> criteria_vector = get_criteria_vector(alternatives);
	std::vector<std::vector<std::vector<double>>> compare_vector;
	for (size_t i = 0; i < 5; ++i)
	{
		std::vector<std::vector<double>> b;
		double number;
		for (size_t j = 0; j < 4; ++j)
		{
			std::vector<double> c;
			for (size_t k = 0; k < 4; ++k)
			{
				if (i != 4) number = (get_mark(criteria_vector[i][j] - criteria_vector[i][k])); // i - criteria name; j, k - alternatives name 
				else number = (get_mark(vec_weight_criteria[j] - vec_weight_criteria[k]));
				c.push_back(number);
			}
			b.push_back(c);
			c.clear();
		}
		compare_vector.push_back(b);
		b.clear();
	}
	std::vector<std::vector<double>> normalized_by_string_vectors;
	for (size_t i = 0; i < 5; ++i)
		normalized_by_string_vectors.push_back(normalise_matrix_by_string(compare_vector[i]));
	for (int i = 0; i < 5; ++i)
	{
		if (i != 4) std::cout << "Compare matrix for " << criteria[i].name << ":\n";
		else std::cout << "Matrix for assesment criteria priorities\n";
		string_delimetr(98, '_');
		for (size_t j = 0; j < alternatives.size(); ++j)
		{
			if (i != 4) std::cout << "	     " << alternatives[j].name;
			else std::cout << "	     " << criteria[j].name;
		}
		std::cout << "	     Normalized by string";
		line;
		for (int j = 0; j < 4; ++j)
		{
			if (i != 4) std::cout << alternatives[j].name;
			else std::cout << criteria[j].name;
			for (int k = 0; k < 5; ++k)
			{
				if (i != 4)
				{
					if (k != 4) std::cout << "	     " << compare_vector[i][j][k];
					else std::cout << std::fixed << std::setprecision(3) << "		     " << normalized_by_string_vectors[i][j];
				}
				else
				{
					if (k != 4) std::cout << "	" << compare_vector[i][j][k] << "	";
					else std::cout << std::fixed << std::setprecision(3) << "	    " << normalized_by_string_vectors[i][j];
				}
			}
			std::cout << std::fixed << std::setprecision(2);
			line;
		}
		string_delimetr(98, '_');
	}
	std::vector<std::vector<double>> matrix;
	for (int i = 0; i < 4; ++i)
	{
		std::vector<double> b;
		for (int j = 0; j < 4; ++j)
			b.push_back(normalized_by_string_vectors[j][i]);
		matrix.push_back(b);
		b.clear();
	}
	std::vector<double> result_vector = multiply_matrix(matrix, normalized_by_string_vectors[4]);
	line;
	std::cout << std::fixed << std::setprecision(3);
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
			if (j == 0) std::cout << "	|" << matrix[i][j];
			else std::cout << "	" << matrix[i][j];
		if (i == 1) std::cout << "|  x  |" << normalized_by_string_vectors[4][i] << "|  =  |" << result_vector[i] << "|\n";
		else std::cout << "|     |" << normalized_by_string_vectors[4][i] << "|  =  |" << result_vector[i] << "|\n";
	}

	std::cout << "\nReceived vector for alternatives:\n";
	string_delimetr(71, '_');
	for (size_t i = 0; i < alternatives.size(); ++i)
		std::cout << "	     " << alternatives[i].name;
	line;
	std::cout << "Assesment";
	for (size_t i = 0; i < 4; ++i)
		std::cout << "     " << result_vector[i] << "	 ";
	line;
	string_delimetr(71, '_');
	int k = 0;
	double max = 0;
	for (size_t i = 0; i < 4; ++i)
		if (result_vector[i] > max)
		{
			max = result_vector[i];
			k = i;
		}
	std::cout << "\n|> Alternative with the best assesment - " << alternatives[k].name << " (Assesment: " << max << ")";
	return k;
}


int main()
{
	std::vector<alternative> alternatives;
	std::vector<criterion> criteria;
	std::vector<double> vec_weight_criteria = {6, 2, 4, 8};
	init1(alternatives);
	init2(criteria);
	try
	{
		print_vec_weight_criteria(criteria, vec_weight_criteria);
	}
	catch (const std::exception& ex)
	{
		std::cout << ex.what();
	}
	int report[4];
	report[0] = method_1(alternatives, criteria, vec_weight_criteria);
	report[1] = method_2(alternatives, criteria);
	report[2] = method_3(alternatives, criteria);
	report[3] = method_4(alternatives, criteria, vec_weight_criteria);
	std::cout << "\n\n\n\t\t\t\tFINAL ANSWERS:\n\n";
	std::cout << "	Method						     Answer\n";
	string_delimetr(71, '-');
	std::cout << "\t1) Method for replacing criteria by limitations:" << std::setw(12) << std::right << alternatives[report[0]].name << std::endl;
	std::cout << "\t2) The formation and narrowing of the Pareto set:" << std::setw(10) << std::right << alternatives[report[1]].name << std::endl;
	std::cout << "\t3) Weighting and combining criteria:" << std::setw(23) << std::right << alternatives[report[2]].name << std::endl;
	std::cout << "\t4) Hierarchy Analysis Method:" << std::setw(30) << std::right << alternatives[report[3]].name << std::endl;
}