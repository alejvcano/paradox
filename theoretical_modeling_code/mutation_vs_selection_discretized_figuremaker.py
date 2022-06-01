# This script specifies the mutation rates, selection coefficients, and nominal probability distribution for the sake of
# calculating the distributions of de novo and fixed mutations; plots are made and saved for all three distributions,
# including regressions, eigenvectors (lengths=standard deviations along 1st and 2nd principal components), and Pearson
# correlations:

import numpy as np
import matplotlib.pyplot as plt
import os

project_parent_directory = '/DIRECTORY_NAME/'

def mutation_selection_parameters():
    # Define range of relative mutation rates and selection coefficients:
    min_mutation = 1
    max_mutation = 5
    min_selection = min_mutation
    max_selection = max_mutation
    values_per_axis = 3

    # Assign relative mutation rates to one axis in 2D plot:
    number_of_mutation_intervals = np.subtract(values_per_axis, 1)
    mutation_range = np.subtract(max_mutation, min_mutation)
    mutation_interval = np.divide(mutation_range, number_of_mutation_intervals)
    half_mut_interval = np.divide(mutation_interval, 2)
    m_ticks = np.arange(min_mutation, max_mutation+mutation_interval, mutation_interval)
    mutation = np.array(m_ticks)

    # Assign relative selection coefficients to one axis in 2D plot:
    number_of_selection_intervals = np.subtract(values_per_axis, 1)
    selection_range = np.subtract(max_selection, min_selection)
    selection_interval = np.divide(selection_range, number_of_selection_intervals)
    half_sel_interval = np.divide(selection_interval, 2)
    s_ticks = np.arange(min_selection, max_selection+selection_interval, selection_interval)
    selection_coefficient_array = np.array(s_ticks)
    selection = np.array([selection_coefficient_array]).T

    # Compute relative rate of fixation (product of mutation and selection):
    ms = np.multiply(mutation, selection)

    # Define ranges of axes in 2D plot:
    low_limit_mut = np.subtract(min_mutation, half_mut_interval)
    high_limit_mut = np.add(max_mutation, half_mut_interval)
    low_limit_sel = np.subtract(min_selection, half_sel_interval)
    high_limit_sel = np.add(max_selection, half_sel_interval)
    m_lim = np.array([low_limit_mut, high_limit_mut])
    s_lim = np.array([low_limit_sel, high_limit_sel])

    # Define the nominal distribution by specifying the number of mutation-rate/selection-coefficient 
    # joint values populated by a mutation class, and the proportion of mutations in that class relative 
    # to total mutations:
    high = np.array([1, 0, 1])
    med = np.array([1, 1, 0])
    low = np.array([0, 1, 1])
    nominal_binary = np.array([low, med, high])
    how_many_mutation_classes = sum(sum(nominal_binary))
    N = np.divide(1, how_many_mutation_classes)
    nominal_proportions = np.multiply(N, nominal_binary)
    
    # Compute axis tick-labels using the range and intervals of mutation and selection defined above:
    mut_list = []
    sel_list = []
    v = 1
    mut_list.append(r'1')
    sel_list.append("1")
    while len(mut_list) < values_per_axis:
        v += int(mutation_interval)
        string_v = str(v)
        next_label_m = r''+string_v+''
        next_label_s = ""+string_v+""
        mut_list.append(next_label_m)
        sel_list.append(next_label_s)
    labels_m = tuple(mut_list)
    labels_s = tuple(sel_list)
    print(labels_m, labels_s)

    return mutation, mutation_interval, m_ticks, selection, selection_interval, s_ticks, selection_coefficient_array, \
           ms, low_limit_mut, high_limit_mut, low_limit_sel, high_limit_sel, m_lim, s_lim, \
           labels_m, labels_s, values_per_axis, nominal_proportions

# Parameter values from the above function, to be used for making figure plots and statistics (PCA, regression lines):
parameters = mutation_selection_parameters()
m = parameters[0]
mutation_interval = parameters[1]
m_ticks = parameters[2]
s = parameters[3]
selection_interval = parameters[4]
s_ticks = parameters[5]
selection_values = parameters[6]
ms = parameters[7]
low_limit_mut = parameters[8]
high_limit_mut = parameters[9]
low_limit_sel = parameters[10]
high_limit_sel = parameters[11]
m_lim = parameters[12]
s_lim = parameters[13]
labels_m = parameters[14]
labels_s = parameters[15]
values_per_axis = parameters[16]
nominal_dist = parameters[17]

# Make folder in which output files will be sorted and saved:
project_output_folder = project_parent_directory + f"/{values_per_axis}_by_{values_per_axis}/"
if not os.path.exists(project_output_folder):
       os.makedirs(project_output_folder)

def mutation_selection_figmaker_nominal():
       # Define and plot the nominal distribution, with relative mutation rate on horizontal axis and
       # relative selection coefficient on vertical axis:
       resized_nom = np.multiply(nominal_dist, 16000)
       sizes = np.array(resized_nom)
       mutation, selection = np.meshgrid(m, s)
       fig_nom, ax_nom = plt.subplots(figsize=(5,5))
       ax_nom.scatter(mutation, selection, s=sizes, c='0.75')
       ax_nom.set(xlim=m_lim, ylim=s_lim)
       plt.subplots_adjust(bottom=0.12, left=0.12)
       plt.xticks(ticks=m_ticks, labels=(labels_m), fontsize=28)
       plt.yticks(ticks=s_ticks, labels=(labels_s), rotation='vertical', verticalalignment='center', fontsize=28)

       # Calculate prevalence-weighted average mutation rate
       total_prevalence = sum(sum(nominal_dist))
       m_prevalences = np.multiply(m, nominal_dist)
       total_m = sum(sum(m_prevalences))
       m_weighted_ave = np.divide(total_m, total_prevalence)

       # Calculate prevalence-weighted average selection coefficients
       transposed_dist = (nominal_dist).T
       nominal_s_dist = sum(transposed_dist)
       nominal_s_dist_trans = (nominal_s_dist).T
       s_prevalences = np.multiply(selection_values, nominal_s_dist_trans)
       total_s = sum(s_prevalences)
       s_weighted_ave = np.divide(total_s, total_prevalence)

       # Calculate covariance of mutation and selection
       ms_prevalences = np.multiply(ms, nominal_dist)
       summed_ms_prev = sum(sum(ms_prevalences))
       mean_of_ms = np.divide(summed_ms_prev, total_prevalence)
       product_of_means = np.multiply(m_weighted_ave, s_weighted_ave)
       cov_ms_raw = np.subtract(mean_of_ms, product_of_means)
       if np.abs(cov_ms_raw) < 0.00000001:
              cov_ms = 0.00
       else:
              cov_ms = round(cov_ms_raw, 10)

       # Calculate variance of mutation rate
       m_squared = np.multiply(m, m)
       m_squared_prev = np.multiply(m_squared, nominal_dist)
       m_squared_weight = np.divide(sum(sum(m_squared_prev)), total_prevalence)
       var_m = np.subtract(m_squared_weight, np.multiply(m_weighted_ave, m_weighted_ave))

       # Calculate variance of selection coefficient
       s_squared_prev = np.multiply((np.multiply(selection_values, selection_values)), nominal_s_dist_trans)
       s_squared_weight = np.divide(sum(s_squared_prev), total_prevalence)
       var_s = np.subtract(s_squared_weight, np.multiply(s_weighted_ave, s_weighted_ave))

       # Calculate and compare eigenvalues for principal component analysis:
       neg_b = np.add(var_m, var_s)
       b_squared = np.multiply(neg_b, neg_b)
       product_vars = np.multiply(var_m, var_s)
       cov_squared = np.multiply(cov_ms, cov_ms)
       c = np.subtract(product_vars, cov_squared)
       c4 = np.multiply(4, c)
       b_squared_minus_4ac = np.subtract(b_squared, c4)
       square_root_quad_form = np.sqrt(b_squared_minus_4ac)
       numerator_plus = np.add(neg_b, square_root_quad_form)
       numerator_minus = np.subtract(neg_b, square_root_quad_form)
       eigenvalue_plus = np.divide(numerator_plus, 2)
       eigenvalue_minus = np.divide(numerator_minus, 2)
       if eigenvalue_plus < 0:
              abs_eigenvalue_plus = np.subtract(0, eigenvalue_plus)
       else:
              abs_eigenvalue_plus = eigenvalue_plus
       if eigenvalue_minus < 0:
              abs_eigenvalue_minus = np.subtract(0, eigenvalue_minus)
       else:
              abs_eigenvalue_minus = eigenvalue_minus
       if abs_eigenvalue_plus >= abs_eigenvalue_minus:
              first_eigenvalue = abs_eigenvalue_plus
              second_eigenvalue = abs_eigenvalue_minus
       else:
              first_eigenvalue = abs_eigenvalue_minus
              second_eigenvalue = abs_eigenvalue_plus

       # Calculate size of principal component axes (with lengths = their standard deviations):
       diff_variances_pc1 = np.subtract(first_eigenvalue, var_m)
       stdev_pc1 = np.sqrt(first_eigenvalue)
       stdev_pc2 = np.sqrt(second_eigenvalue)
       half_stdev_pc1 = np.divide(stdev_pc1, 2)
       half_stdev_pc2 = np.divide(stdev_pc2, 2)

       # Calculate slopes and positions of principal component axes and plot them:
       if cov_ms == 0:
              if round(var_m, 10) >= round(var_s, 10):
                     slope_pc1 = 0
                     slope_pc2 = np.inf
                     plt.hlines(y=s_weighted_ave, xmin=m_weighted_ave-half_stdev_pc1,
                                xmax=m_weighted_ave+half_stdev_pc1, color='k', linewidth=3, label="PC1 regression")
                     plt.vlines(x=m_weighted_ave, ymin=s_weighted_ave-half_stdev_pc2,
                                ymax=s_weighted_ave+half_stdev_pc2, color='k', linewidth=3, label="PC1 regression")
              else:
                     slope_pc1 = np.inf
                     slope_pc2 = 0
                     plt.hlines(y=s_weighted_ave, xmin=m_weighted_ave-half_stdev_pc2,
                                xmax=m_weighted_ave+half_stdev_pc2, color='k', linewidth=3, label="PC1 regression")
                     plt.vlines(x=m_weighted_ave, ymin=s_weighted_ave-half_stdev_pc1,
                                ymax=s_weighted_ave+half_stdev_pc1, color='k', linewidth=3, label="PC1 regression")
       else:
              slope_pc1 = np.divide(diff_variances_pc1, cov_ms)
              slope_pc2 = np.divide(-1, slope_pc1)
              mx_pc1 = np.multiply(slope_pc1, m_weighted_ave)
              mx_pc2 = np.multiply(slope_pc2, m_weighted_ave)
              intercept_pc1 = np.subtract(s_weighted_ave, mx_pc1)
              intercept_pc2 = np.subtract(s_weighted_ave, mx_pc2)
              intercept_minus_mean_pc1 = np.subtract(intercept_pc1, s_weighted_ave)
              intercept_minus_mean_pc2 = np.subtract(intercept_pc2, s_weighted_ave)
              if intercept_minus_mean_pc1 < 0:
                     intercept_mean_vertical_distance_pc1 = np.subtract(0, intercept_minus_mean_pc1)
              else:
                     intercept_mean_vertical_distance_pc1 = intercept_minus_mean_pc1
              if intercept_minus_mean_pc2 < 0:
                     intercept_mean_vertical_distance_pc2 = np.subtract(0, intercept_minus_mean_pc2)
              else:
                     intercept_mean_vertical_distance_pc2 = intercept_minus_mean_pc2
              opp_over_adj_pc1 = np.divide(intercept_mean_vertical_distance_pc1, m_weighted_ave)
              opp_over_adj_pc2 = np.divide(intercept_mean_vertical_distance_pc2, m_weighted_ave)
              angle_pc1 = np.arctan(opp_over_adj_pc1)
              angle_pc2 = np.arctan(opp_over_adj_pc2)
              cos_angle_pc1 = np.cos(angle_pc1)
              cos_angle_pc2 = np.cos(angle_pc2)
              horizontal_dist_half_stdev_pc1 = np.multiply(cos_angle_pc1, half_stdev_pc1)
              horizontal_dist_half_stdev_pc2 = np.multiply(cos_angle_pc2, half_stdev_pc2)
              m_low_pc1 = np.subtract(m_weighted_ave, horizontal_dist_half_stdev_pc1)
              m_high_pc1 = np.add(m_weighted_ave, horizontal_dist_half_stdev_pc1)
              m_low_pc2 = np.subtract(m_weighted_ave, horizontal_dist_half_stdev_pc2)
              m_high_pc2 = np.add(m_weighted_ave, horizontal_dist_half_stdev_pc2)
              m_range_pc1 = np.linspace(m_low_pc1, m_high_pc1, 100)
              m_range_pc2 = np.linspace(m_low_pc2, m_high_pc2, 100)
              regression_pc1 = slope_pc1*m_range_pc1+intercept_pc1
              regression_pc2 = slope_pc2*m_range_pc2+intercept_pc2
              plt.plot(m_range_pc1, regression_pc1, '-k', linewidth=3, label="PC1 regression")
              plt.plot(m_range_pc2, regression_pc2, '-k', linewidth=3, label="PC2 regression")

       # Linear regressions (slope & intercept for selection=f(mutation) and mutation=f(selection)):
       slope_function_of_mutation = np.divide(cov_ms, var_m)
       slope_function_of_selection = np.divide(cov_ms, var_s)
       mx_function_of_mutation = np.multiply(slope_function_of_mutation, m_weighted_ave)
       my_function_of_selection = np.multiply(slope_function_of_selection, s_weighted_ave)
       intercept_function_of_mutation = np.subtract(s_weighted_ave, mx_function_of_mutation)
       intercept_function_of_selection = np.subtract(m_weighted_ave, my_function_of_selection)

       # Plot selection as a function of mutation:
       plt.plot(m_lim, slope_function_of_mutation * m_lim + intercept_function_of_mutation, '-k', linewidth=1.5,
                label="Mutation regression")

       # Plot mutation as a function of selection:
       if slope_function_of_selection == 0:
              plt.vlines(x=m_weighted_ave, ymin=low_limit_sel, ymax=high_limit_sel, colors='k', linestyles='dotted',
                         linewidth=1.5, label="Selection regression")
       else:
              plt.plot(m_lim, (m_lim - intercept_function_of_selection) / slope_function_of_selection, ':k',
                       linewidth=1.5, label="Selection regression")

       # Calculate correlation coefficient
       prod_of_vars = np.multiply(var_m, var_s)
       sqrt_of_vars = np.sqrt(prod_of_vars)
       rho_ms = np.divide(cov_ms, sqrt_of_vars)
       if np.abs(rho_ms) < 0.00000001:
              r_nom = 0
              rho_ms_text = "r = 0"
       else:
              r_nom = "{:.4f}".format(rho_ms)
              rho_ms_text = "r = {:.2f}".format(rho_ms)
       print("nominal correlation ", rho_ms_text)

       # Position and show correlation coefficient on figure graph:
       half_height = np.divide(high_limit_mut, 2)
       if values_per_axis == 2:
              dist_from_top = np.divide(selection_interval, 4.4)
              dist_from_bottom = np.divide(selection_interval, 12)
              dist_from_left_side = np.divide(mutation_interval, 12)
              if np.abs(rho_ms) < 0.00000001:
                     dist_from_right_side = np.divide(mutation_interval, 1.7)
              else:
                     dist_from_right_side = np.divide(mutation_interval, 1.05)
       else:
              dist_from_top = np.divide(selection_interval, 3)
              dist_from_bottom = np.divide(selection_interval, 8)
              dist_from_left_side = np.divide(mutation_interval, 8)
              if np.abs(rho_ms) < 0.00000001:
                     dist_from_right_side = np.divide(mutation_interval, 1.2)
              else:
                     dist_from_right_side = np.divide(mutation_interval, 0.7)
       if s_weighted_ave <= half_height and slope_pc1 <= 0:
              rho_x_coordinate = np.subtract(high_limit_mut, dist_from_right_side)
              rho_y_coordinate = np.subtract(high_limit_sel, dist_from_top)
       elif s_weighted_ave <= half_height and slope_pc1 > 0:
              rho_x_coordinate = np.add(low_limit_mut, dist_from_left_side)
              rho_y_coordinate = np.subtract(high_limit_sel, dist_from_top)
       elif s_weighted_ave > half_height and slope_pc1 <= 0:
              rho_x_coordinate = np.add(low_limit_mut, dist_from_left_side)
              rho_y_coordinate = np.add(low_limit_sel, dist_from_bottom)
       else:
              rho_x_coordinate = np.subtract(high_limit_mut, dist_from_right_side)
              rho_y_coordinate = np.add(low_limit_sel, dist_from_bottom)
       if np.abs(rho_ms) < 0.00000001:
              plt.text(rho_x_coordinate, rho_y_coordinate, "r = 0", fontdict=None, fontsize=28)
       else:
              plt.text(rho_x_coordinate, rho_y_coordinate, rho_ms_text, fontdict=None, fontsize=28)

       # Save figure files to destination folder:
       nominal_outputs = project_output_folder + '/nominal_plots/'
       if not os.path.exists(nominal_outputs):
              os.makedirs(nominal_outputs)
       n = 1
       while os.path.exists(nominal_outputs+f"fig_nom_{n}.png"):
              n += 1
       fig_nom_file = plt.savefig(nominal_outputs+"fig_nom_{}.png".format(n))
       plt.close(fig_nom_file)

       # Show the plot (de-comment the following line if you want the script to open the plot when it runs):
       # plt.show()

       return fig_nom_file, r_nom

def mutation_selection_figmaker_de_novo():
       # Calculate and plot the de novo distribution, with relative mutation rate on horizontal axis and
       # relative selection coefficient on vertical axis:
       m_weighted = np.multiply(m, nominal_dist)
       total_area = np.sum(m_weighted)
       de_novo_proportions = np.divide(m_weighted, total_area)
       resized_de_novo = np.multiply(de_novo_proportions, 16000)
       adjusted_sizes = np.array(resized_de_novo)
       mutation, selection = np.meshgrid(m, s)
       fig_de_novo, ax_de_novo = plt.subplots(figsize=(5,5))
       ax_de_novo.scatter(mutation, selection, s=adjusted_sizes, c='0.75')
       ax_de_novo.set(xlim=m_lim, ylim=s_lim)
       plt.subplots_adjust(bottom=0.12, left=0.12)
       plt.xticks(ticks=m_ticks, labels=(labels_m), fontsize=28)
       plt.yticks(ticks=s_ticks, labels=(labels_s), rotation='vertical', verticalalignment='center', fontsize=28)

       # Calculate prevalence-weighted average mutation rate
       total_prevalence = sum(sum(de_novo_proportions))
       m_prevalences = np.multiply(mutation, de_novo_proportions)
       total_m = sum(sum(m_prevalences))
       m_weighted_ave = np.divide(total_m, total_prevalence)

       # Calculate prevalence-weighted average selection coefficients
       transposed_dist = (de_novo_proportions).T
       de_novo_s_dist = sum(transposed_dist)
       de_novo_s_dist_trans = (de_novo_s_dist).T
       s_prevalences = np.multiply(selection_values, de_novo_s_dist_trans)
       total_s = sum(s_prevalences)
       s_weighted_ave = np.divide(total_s, total_prevalence)

       # Calculate covariance of mutation and selection
       ms_prevalences = np.multiply(ms, de_novo_proportions)
       summed_ms_prev = sum(sum(ms_prevalences))
       mean_of_ms = np.divide(summed_ms_prev, total_prevalence)
       product_of_means = np.multiply(m_weighted_ave, s_weighted_ave)
       cov_ms_raw = np.subtract(mean_of_ms, product_of_means)
       if np.abs(cov_ms_raw) < 0.00000001:
              cov_ms = 0.00
       else:
              cov_ms = round(cov_ms_raw, 10)

       # Calculate variance of mutation rate
       m_squared = np.multiply(mutation, mutation)
       m_squared_prev = np.multiply(m_squared, de_novo_proportions)
       m_squared_weight = np.divide(sum(sum(m_squared_prev)), total_prevalence)
       var_m = np.subtract(m_squared_weight, np.multiply(m_weighted_ave, m_weighted_ave))

       # Calculate variance of selection coefficient
       s_squared_prev = np.multiply((np.multiply(selection_values, selection_values)), de_novo_s_dist_trans)
       s_squared_weight = np.divide(sum(s_squared_prev), total_prevalence)
       var_s = np.subtract(s_squared_weight, np.multiply(s_weighted_ave, s_weighted_ave))

       # Calculate and compare eigenvalues for principal component analysis:
       neg_b = np.add(var_m, var_s)
       b_squared = np.multiply(neg_b, neg_b)
       product_vars = np.multiply(var_m, var_s)
       cov_squared = np.multiply(cov_ms, cov_ms)
       c = np.subtract(product_vars, cov_squared)
       c4 = np.multiply(4, c)
       b_squared_minus_4ac = np.subtract(b_squared, c4)
       square_root_quad_form = np.sqrt(b_squared_minus_4ac)
       numerator_plus = np.add(neg_b, square_root_quad_form)
       numerator_minus = np.subtract(neg_b, square_root_quad_form)
       eigenvalue_plus = np.divide(numerator_plus, 2)
       eigenvalue_minus = np.divide(numerator_minus, 2)
       if eigenvalue_plus < 0:
              abs_eigenvalue_plus = np.subtract(0, eigenvalue_plus)
       else:
              abs_eigenvalue_plus = eigenvalue_plus
       if eigenvalue_minus < 0:
              abs_eigenvalue_minus = np.subtract(0, eigenvalue_minus)
       else:
              abs_eigenvalue_minus = eigenvalue_minus
       if abs_eigenvalue_plus >= abs_eigenvalue_minus:
              first_eigenvalue = abs_eigenvalue_plus
              second_eigenvalue = abs_eigenvalue_minus
       else:
              first_eigenvalue = abs_eigenvalue_minus
              second_eigenvalue = abs_eigenvalue_plus

       # Calculate size of principal component axes (with lengths = their standard deviations):
       diff_variances_pc1 = np.subtract(first_eigenvalue, var_m)
       stdev_pc1 = np.sqrt(first_eigenvalue)
       stdev_pc2 = np.sqrt(second_eigenvalue)
       half_stdev_pc1 = np.divide(stdev_pc1, 2)
       half_stdev_pc2 = np.divide(stdev_pc2, 2)

       # Calculate slopes and positions of principal component axes and plot them:
       if cov_ms == 0:
              if round(var_m, 10) >= round(var_s, 10):
                     slope_pc1 = 0
                     slope_pc2 = np.inf
                     plt.hlines(y=s_weighted_ave, xmin=m_weighted_ave - half_stdev_pc1,
                                xmax=m_weighted_ave + half_stdev_pc1, color='k', linewidth=3, label="PC1 regression")
                     plt.vlines(x=m_weighted_ave, ymin=s_weighted_ave - half_stdev_pc2,
                                ymax=s_weighted_ave + half_stdev_pc2, color='k', linewidth=3, label="PC1 regression")
              else:
                     slope_pc1 = np.inf
                     slope_pc2 = 0
                     plt.hlines(y=s_weighted_ave, xmin=m_weighted_ave - half_stdev_pc2,
                                xmax=m_weighted_ave + half_stdev_pc2, color='k', linewidth=3, label="PC1 regression")
                     plt.vlines(x=m_weighted_ave, ymin=s_weighted_ave - half_stdev_pc1,
                                ymax=s_weighted_ave + half_stdev_pc1, color='k', linewidth=3, label="PC1 regression")
       else:
              slope_pc1 = np.divide(diff_variances_pc1, cov_ms)
              slope_pc2 = np.divide(-1, slope_pc1)
              mx_pc1 = np.multiply(slope_pc1, m_weighted_ave)
              mx_pc2 = np.multiply(slope_pc2, m_weighted_ave)
              intercept_pc1 = np.subtract(s_weighted_ave, mx_pc1)
              intercept_pc2 = np.subtract(s_weighted_ave, mx_pc2)
              intercept_minus_mean_pc1 = np.subtract(intercept_pc1, s_weighted_ave)
              intercept_minus_mean_pc2 = np.subtract(intercept_pc2, s_weighted_ave)
              if intercept_minus_mean_pc1 < 0:
                     intercept_mean_vertical_distance_pc1 = np.subtract(0, intercept_minus_mean_pc1)
              else:
                     intercept_mean_vertical_distance_pc1 = intercept_minus_mean_pc1
              if intercept_minus_mean_pc2 < 0:
                     intercept_mean_vertical_distance_pc2 = np.subtract(0, intercept_minus_mean_pc2)
              else:
                     intercept_mean_vertical_distance_pc2 = intercept_minus_mean_pc2
              opp_over_adj_pc1 = np.divide(intercept_mean_vertical_distance_pc1, m_weighted_ave)
              opp_over_adj_pc2 = np.divide(intercept_mean_vertical_distance_pc2, m_weighted_ave)
              angle_pc1 = np.arctan(opp_over_adj_pc1)
              angle_pc2 = np.arctan(opp_over_adj_pc2)
              cos_angle_pc1 = np.cos(angle_pc1)
              cos_angle_pc2 = np.cos(angle_pc2)
              horizontal_dist_half_stdev_pc1 = np.multiply(cos_angle_pc1, half_stdev_pc1)
              horizontal_dist_half_stdev_pc2 = np.multiply(cos_angle_pc2, half_stdev_pc2)
              m_low_pc1 = np.subtract(m_weighted_ave, horizontal_dist_half_stdev_pc1)
              m_high_pc1 = np.add(m_weighted_ave, horizontal_dist_half_stdev_pc1)
              m_low_pc2 = np.subtract(m_weighted_ave, horizontal_dist_half_stdev_pc2)
              m_high_pc2 = np.add(m_weighted_ave, horizontal_dist_half_stdev_pc2)
              m_range_pc1 = np.linspace(m_low_pc1, m_high_pc1, 100)
              m_range_pc2 = np.linspace(m_low_pc2, m_high_pc2, 100)
              regression_pc1 = slope_pc1*m_range_pc1+intercept_pc1
              regression_pc2 = slope_pc2*m_range_pc2+intercept_pc2
              plt.plot(m_range_pc1, regression_pc1, '-k', linewidth=3, label="PC1 regression")
              plt.plot(m_range_pc2, regression_pc2, '-k', linewidth=3, label="PC2 regression")

       # Linear regressions (slope & intercept for selection=f(mutation) and mutation=f(selection)):
       slope_function_of_mutation = np.divide(cov_ms, var_m)
       slope_function_of_selection = np.divide(cov_ms, var_s)
       mx_function_of_mutation = np.multiply(slope_function_of_mutation, m_weighted_ave)
       my_function_of_selection = np.multiply(slope_function_of_selection, s_weighted_ave)
       intercept_function_of_mutation = np.subtract(s_weighted_ave, mx_function_of_mutation)
       intercept_function_of_selection = np.subtract(m_weighted_ave, my_function_of_selection)

       # Plot selection as a function of mutation:
       plt.plot(m_lim, slope_function_of_mutation * m_lim + intercept_function_of_mutation, '-k', linewidth=1.5,
                label="Mutation regression")

       # Plot mutation as a function of selection:
       if slope_function_of_selection == 0:
              plt.vlines(x=m_weighted_ave, ymin=low_limit_sel, ymax=high_limit_sel, colors='k', linestyles='dotted',
                         linewidth=1.5, label="Selection regression")
       else:
              plt.plot(m_lim, (m_lim - intercept_function_of_selection) / slope_function_of_selection, ':k',
                       linewidth=1.5, label="Selection regression")

       # Calculate correlation coefficient
       prod_of_vars = np.multiply(var_m, var_s)
       sqrt_of_vars = np.sqrt(prod_of_vars)
       rho_ms = np.divide(cov_ms, sqrt_of_vars)
       if np.abs(rho_ms) < 0.00000001:
              r_dn = 0
              rho_ms_text = "r = 0"
       else:
              r_dn = "{:.4f}".format(rho_ms)
              rho_ms_text = "r = {:.2f}".format(rho_ms)
       print("de novo correlation ", rho_ms_text)

       # Position and show correlation coefficient on figure graph:
       half_height = np.divide(high_limit_mut, 2)
       if values_per_axis == 2:
              dist_from_top = np.divide(selection_interval, 4.4)
              dist_from_bottom = np.divide(selection_interval, 12)
              dist_from_left_side = np.divide(mutation_interval, 12)
              if np.abs(rho_ms) < 0.00000001:
                     dist_from_right_side = np.divide(mutation_interval, 1.7)
              else:
                     dist_from_right_side = np.divide(mutation_interval, 1.05)
       else:
              dist_from_top = np.divide(selection_interval, 3)
              dist_from_bottom = np.divide(selection_interval, 8)
              dist_from_left_side = np.divide(mutation_interval, 8)
              if np.abs(rho_ms) < 0.00000001:
                     dist_from_right_side = np.divide(mutation_interval, 1.2)
              else:
                     dist_from_right_side = np.divide(mutation_interval, 0.7)
       if s_weighted_ave <= half_height and slope_pc1 <= 0:
              rho_x_coordinate = np.subtract(high_limit_mut, dist_from_right_side)
              rho_y_coordinate = np.subtract(high_limit_sel, dist_from_top)
       elif s_weighted_ave <= half_height and slope_pc1 > 0:
              rho_x_coordinate = np.add(low_limit_mut, dist_from_left_side)
              rho_y_coordinate = np.subtract(high_limit_sel, dist_from_top)
       elif s_weighted_ave > half_height and slope_pc1 <= 0:
              rho_x_coordinate = np.add(low_limit_mut, dist_from_left_side)
              rho_y_coordinate = np.add(low_limit_sel, dist_from_bottom)
       else:
              rho_x_coordinate = np.subtract(high_limit_mut, dist_from_right_side)
              rho_y_coordinate = np.add(low_limit_sel, dist_from_bottom)
       if np.abs(rho_ms) < 0.00000001:
              plt.text(rho_x_coordinate, rho_y_coordinate, "r = 0", fontdict=None, fontsize=28)
       else:
              plt.text(rho_x_coordinate, rho_y_coordinate, rho_ms_text, fontdict=None, fontsize=28)

       # Save figure files to destination folder:
       de_novo_outputs = project_output_folder + '/de_novo_plots/'
       if not os.path.exists(de_novo_outputs):
              os.makedirs(de_novo_outputs)
       n = 1
       while os.path.exists(de_novo_outputs+f"fig_de_novo_{n}.png"):
              n += 1
       fig_de_novo_file = plt.savefig(de_novo_outputs+"fig_de_novo_{}.png".format(n))
       plt.close(fig_de_novo_file)

       # Show the plot (de-comment the following line if you want the script to open the plot when it runs):
       # plt.show()

       return fig_de_novo_file, r_dn

def mutation_selection_figmaker_fixed():
       # Calculate and plot the fixed distribution, with relative mutation rate on horizontal axis and
       # relative selection coefficient on vertical axis:
       ms_weighted = np.multiply(ms, nominal_dist)
       total_area = np.sum(ms_weighted)
       fixed_proportions = np.divide(ms_weighted, total_area)
       resized_fixed = np.multiply(fixed_proportions, 16000)
       adjusted_sizes = np.array(resized_fixed)
       mutation, selection = np.meshgrid(m, s)
       fig_fixed, ax_fixed = plt.subplots(figsize=(5,5))
       ax_fixed.scatter(mutation, selection, s=adjusted_sizes, c='0.75')
       ax_fixed.set(xlim=m_lim, ylim=s_lim)
       plt.subplots_adjust(bottom=0.12, left=0.12)
       plt.xticks(ticks=m_ticks, labels=(labels_m), fontsize=28)
       plt.yticks(ticks=s_ticks, labels=(labels_s), rotation='vertical', verticalalignment='center', fontsize=28)

       # Calculate prevalence-weighted average mutation rate
       total_prevalence = sum(sum(fixed_proportions))
       m_prevalences = np.multiply(mutation, fixed_proportions)
       total_m = sum(sum(m_prevalences))
       m_weighted_ave = np.divide(total_m, total_prevalence)

       # Calculate prevalence-weighted average selection coefficients
       transposed_dist = (fixed_proportions).T
       fixed_s_dist = sum(transposed_dist)
       fixed_s_dist_trans = (fixed_s_dist).T
       s_prevalences = np.multiply(selection_values, fixed_s_dist_trans)
       total_s = sum(s_prevalences)
       s_weighted_ave = np.divide(total_s, total_prevalence)

       # Calculate covariance of mutation and selection
       ms_prevalences = np.multiply(ms, fixed_proportions)
       summed_ms_prev = sum(sum(ms_prevalences))
       mean_of_ms = np.divide(summed_ms_prev, total_prevalence)
       product_of_means = np.multiply(m_weighted_ave, s_weighted_ave)
       cov_ms_raw = np.subtract(mean_of_ms, product_of_means)
       if np.abs(cov_ms_raw) < 0.00000001:
              cov_ms = 0.00
       else:
              cov_ms = round(cov_ms_raw, 10)

       # Calculate variance of mutation rate
       m_squared = np.multiply(mutation, mutation)
       m_squared_prev = np.multiply(m_squared, fixed_proportions)
       m_squared_weight = np.divide(sum(sum(m_squared_prev)), total_prevalence)
       var_m = np.subtract(m_squared_weight, np.multiply(m_weighted_ave, m_weighted_ave))

       # Calculate variance of selection coefficient
       s_squared_prev = np.multiply((np.multiply(selection_values, selection_values)), fixed_s_dist_trans)
       s_squared_weight = np.divide(sum(s_squared_prev), total_prevalence)
       var_s = np.subtract(s_squared_weight, np.multiply(s_weighted_ave, s_weighted_ave))

       # Calculate and compare eigenvalues for principal component analysis:
       neg_b = np.add(var_m, var_s)
       b_squared = np.multiply(neg_b, neg_b)
       product_vars = np.multiply(var_m, var_s)
       cov_squared = np.multiply(cov_ms, cov_ms)
       c = np.subtract(product_vars, cov_squared)
       c4 = np.multiply(4, c)
       b_squared_minus_4ac = np.subtract(b_squared, c4)
       square_root_quad_form = np.sqrt(b_squared_minus_4ac)
       numerator_plus = np.add(neg_b, square_root_quad_form)
       numerator_minus = np.subtract(neg_b, square_root_quad_form)
       eigenvalue_plus = np.divide(numerator_plus, 2)
       eigenvalue_minus = np.divide(numerator_minus, 2)
       if eigenvalue_plus < 0:
              abs_eigenvalue_plus = np.subtract(0, eigenvalue_plus)
       else:
              abs_eigenvalue_plus = eigenvalue_plus
       if eigenvalue_minus < 0:
              abs_eigenvalue_minus = np.subtract(0, eigenvalue_minus)
       else:
              abs_eigenvalue_minus = eigenvalue_minus
       if abs_eigenvalue_plus >= abs_eigenvalue_minus:
              first_eigenvalue = abs_eigenvalue_plus
              second_eigenvalue = abs_eigenvalue_minus
       else:
              first_eigenvalue = abs_eigenvalue_minus
              second_eigenvalue = abs_eigenvalue_plus

       # Calculate size of principal component axes (with lengths = their standard deviations):
       diff_variances_pc1 = np.subtract(first_eigenvalue, var_m)
       stdev_pc1 = np.sqrt(first_eigenvalue)
       stdev_pc2 = np.sqrt(second_eigenvalue)
       half_stdev_pc1 = np.divide(stdev_pc1, 2)
       half_stdev_pc2 = np.divide(stdev_pc2, 2)

       # Calculate slopes and positions of principal component axes and plot them:
       if cov_ms == 0:
              if round(var_m, 10) >= round(var_s, 10):
                     slope_pc1 = 0
                     slope_pc2 = np.inf
                     plt.hlines(y=s_weighted_ave, xmin=m_weighted_ave - half_stdev_pc1,
                                xmax=m_weighted_ave + half_stdev_pc1, color='k', linewidth=3, label="PC1 regression")
                     plt.vlines(x=m_weighted_ave, ymin=s_weighted_ave - half_stdev_pc2,
                                ymax=s_weighted_ave + half_stdev_pc2, color='k', linewidth=3, label="PC1 regression")
              else:
                     slope_pc1 = np.inf
                     slope_pc2 = 0
                     plt.hlines(y=s_weighted_ave, xmin=m_weighted_ave - half_stdev_pc2,
                                xmax=m_weighted_ave + half_stdev_pc2, color='k', linewidth=3, label="PC1 regression")
                     plt.vlines(x=m_weighted_ave, ymin=s_weighted_ave - half_stdev_pc1,
                                ymax=s_weighted_ave + half_stdev_pc1, color='k', linewidth=3, label="PC1 regression")
       else:
              slope_pc1 = np.divide(diff_variances_pc1, cov_ms)
              slope_pc2 = np.divide(-1, slope_pc1)
              mx_pc1 = np.multiply(slope_pc1, m_weighted_ave)
              mx_pc2 = np.multiply(slope_pc2, m_weighted_ave)
              intercept_pc1 = np.subtract(s_weighted_ave, mx_pc1)
              intercept_pc2 = np.subtract(s_weighted_ave, mx_pc2)
              intercept_minus_mean_pc1 = np.subtract(intercept_pc1, s_weighted_ave)
              intercept_minus_mean_pc2 = np.subtract(intercept_pc2, s_weighted_ave)
              if intercept_minus_mean_pc1 < 0:
                     intercept_mean_vertical_distance_pc1 = np.subtract(0, intercept_minus_mean_pc1)
              else:
                     intercept_mean_vertical_distance_pc1 = intercept_minus_mean_pc1
              if intercept_minus_mean_pc2 < 0:
                     intercept_mean_vertical_distance_pc2 = np.subtract(0, intercept_minus_mean_pc2)
              else:
                     intercept_mean_vertical_distance_pc2 = intercept_minus_mean_pc2
              opp_over_adj_pc1 = np.divide(intercept_mean_vertical_distance_pc1, m_weighted_ave)
              opp_over_adj_pc2 = np.divide(intercept_mean_vertical_distance_pc2, m_weighted_ave)
              angle_pc1 = np.arctan(opp_over_adj_pc1)
              angle_pc2 = np.arctan(opp_over_adj_pc2)
              cos_angle_pc1 = np.cos(angle_pc1)
              cos_angle_pc2 = np.cos(angle_pc2)
              horizontal_dist_half_stdev_pc1 = np.multiply(cos_angle_pc1, half_stdev_pc1)
              horizontal_dist_half_stdev_pc2 = np.multiply(cos_angle_pc2, half_stdev_pc2)
              m_low_pc1 = np.subtract(m_weighted_ave, horizontal_dist_half_stdev_pc1)
              m_high_pc1 = np.add(m_weighted_ave, horizontal_dist_half_stdev_pc1)
              m_low_pc2 = np.subtract(m_weighted_ave, horizontal_dist_half_stdev_pc2)
              m_high_pc2 = np.add(m_weighted_ave, horizontal_dist_half_stdev_pc2)
              m_range_pc1 = np.linspace(m_low_pc1, m_high_pc1, 100)
              m_range_pc2 = np.linspace(m_low_pc2, m_high_pc2, 100)
              regression_pc1 = slope_pc1*m_range_pc1+intercept_pc1
              regression_pc2 = slope_pc2*m_range_pc2+intercept_pc2
              plt.plot(m_range_pc1, regression_pc1, '-k', linewidth=3, label="PC1 regression")
              plt.plot(m_range_pc2, regression_pc2, '-k', linewidth=3, label="PC2 regression")

       # Linear regressions (slope & intercept for selection=f(mutation) and mutation=f(selection)):
       slope_function_of_mutation = np.divide(cov_ms, var_m)
       slope_function_of_selection = np.divide(cov_ms, var_s)
       mx_function_of_mutation = np.multiply(slope_function_of_mutation, m_weighted_ave)
       my_function_of_selection = np.multiply(slope_function_of_selection, s_weighted_ave)
       intercept_function_of_mutation = np.subtract(s_weighted_ave, mx_function_of_mutation)
       intercept_function_of_selection = np.subtract(m_weighted_ave, my_function_of_selection)

       # Plot selection as a function of mutation:
       plt.plot(m_lim, slope_function_of_mutation * m_lim + intercept_function_of_mutation, '-k', linewidth=1.5,
                label="Mutation regression")

       # Plot mutation as a function of selection:
       if slope_function_of_selection == 0:
              plt.vlines(x=m_weighted_ave, ymin=low_limit_sel, ymax=high_limit_sel, colors='k', linestyles='dotted',
                         linewidth=1.5, label="Selection regression")
       else:
              plt.plot(m_lim, (m_lim - intercept_function_of_selection) / slope_function_of_selection, ':k',
                       linewidth=1.5, label="Selection regression")

       # Calculate correlation coefficient
       prod_of_vars = np.multiply(var_m, var_s)
       sqrt_of_vars = np.sqrt(prod_of_vars)
       rho_ms = np.divide(cov_ms, sqrt_of_vars)
       if np.abs(rho_ms) < 0.00000001:
              r_fix = 0
              rho_ms_text = "r = 0"
       else:
              r_fix = "{:.4f}".format(rho_ms)
              rho_ms_text = "r = {:.2f}".format(rho_ms)
       print("fixed correlation ", rho_ms_text)

       # Position and show correlation coefficient on figure graph:
       half_height = np.divide(high_limit_mut, 2)
       if values_per_axis == 2:
              dist_from_top = np.divide(selection_interval, 4.4)
              dist_from_bottom = np.divide(selection_interval, 12)
              dist_from_left_side = np.divide(mutation_interval, 12)
              if np.abs(rho_ms) < 0.00000001:
                     dist_from_right_side = np.divide(mutation_interval, 1.7)
              else:
                     dist_from_right_side = np.divide(mutation_interval, 1.05)
       else:
              dist_from_top = np.divide(selection_interval, 3)
              dist_from_bottom = np.divide(selection_interval, 8)
              dist_from_left_side = np.divide(mutation_interval, 8)
              if np.abs(rho_ms) < 0.00000001:
                     dist_from_right_side = np.divide(mutation_interval, 1.2)
              else:
                     dist_from_right_side = np.divide(mutation_interval, 0.7)
       if s_weighted_ave <= half_height and slope_pc1 <= 0:
              rho_x_coordinate = np.subtract(high_limit_mut, dist_from_right_side)
              rho_y_coordinate = np.subtract(high_limit_sel, dist_from_top)
       elif s_weighted_ave <= half_height and slope_pc1 > 0:
              rho_x_coordinate = np.add(low_limit_mut, dist_from_left_side)
              rho_y_coordinate = np.subtract(high_limit_sel, dist_from_top)
       elif s_weighted_ave > half_height and slope_pc1 <= 0:
              rho_x_coordinate = np.add(low_limit_mut, dist_from_left_side)
              rho_y_coordinate = np.add(low_limit_sel, dist_from_bottom)
       else:
              rho_x_coordinate = np.subtract(high_limit_mut, dist_from_right_side)
              rho_y_coordinate = np.add(low_limit_sel, dist_from_bottom)
       if np.abs(rho_ms) < 0.00000001:
              plt.text(rho_x_coordinate, rho_y_coordinate, "r = 0", fontdict=None, fontsize=28)
       else:
              plt.text(rho_x_coordinate, rho_y_coordinate, rho_ms_text, fontdict=None, fontsize=28)

       # Save figure files to destination folder:
       fixed_outputs = project_output_folder + '/fixed_plots/'
       if not os.path.exists(fixed_outputs):
              os.makedirs(fixed_outputs)
       n = 1
       while os.path.exists(fixed_outputs+f"fig_fixed_{n}.png"):
              n += 1
       fig_fixed_file = plt.savefig(fixed_outputs+"fig_fixed_{}.png".format(n))
       plt.close(fig_fixed_file)

       # Show the plot (de-comment the following line if you want the script to open the plot when it runs):
       # plt.show()

       return fig_fixed_file, r_fix

if __name__ == '__main__':
    mutation_selection_figmaker_nominal()
    mutation_selection_figmaker_de_novo()
    mutation_selection_figmaker_fixed()
