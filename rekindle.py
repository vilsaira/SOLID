"""
REKINDLE reweighting as published in Tax et al., MRM 2015
Copyright Chantal Tax (taxc@cardiff.ac.uk) and Alexander Leemans
(alexander@isi.uu.nl) .

Python adaptation and SOLID weighting (Sairanen et al., NeuroImage 2019)
implementation by Viljami Sairanen (viljami.sairanen@gmail.com).
"""
import numpy as np
from numpy import dot
import scipy.linalg as linalg

def robust_fit(measurements=None, design_matrix=None, weights=None, 
            iter_max_n_inner_loop=5, iter_max_n_full_loop=20, iter_convergence=1.0e-3,
            kappa_thr_up=6.0, kappa_thr_low=3.0, downweight_outliers=True):
    if measurements is None:
        return -1
    if design_matrix is None:
        return -1
    if weights is None:
        weights = np.ones(measurements.shape)

    # Calculate residuals using REKINDLE approach    
    residuals = _outlier_detection_calculate_residuals(measurements, design_matrix,
                                        iter_max_n_inner_loop, iter_max_n_full_loop, iter_convergence)
    
    # Calculate downweighting factors based on REKINDLE outliers
    C = 1.4826*np.median(np.abs(residuals - np.median(residuals)))
    if C < iter_convergence:
        C = iter_convergence
    outlier_indx = np.abs(residuals) > kappa_thr_up*C
    if downweight_outliers:    
        rekindle_weights = np.abs(residuals)/C
        rekindle_weights[rekindle_weights <= kappa_thr_low] = kappa_thr_low
        rekindle_weights[rekindle_weights > kappa_thr_up] = kappa_thr_up
        rekindle_weights = (rekindle_weights - kappa_thr_low) / (kappa_thr_up - kappa_thr_low)
        rekindle_weights = 1 - rekindle_weights
        weights = weights * rekindle_weights

        beta = iwlls(measurements=measurements, design_matrix=design_matrix,
            weights=weights, iter_max_n=iter_max_n_full_loop, iter_convergence=iter_convergence)
    else:
        # Outliers detected by REKINDLE are excluded instead of downweighted which
        # is the original REKINDLE approach.
        beta = iwlls(measurements=measurements[~outlier_indx], design_matrix=design_matrix[~outlier_indx,:],
            weights=weights[~outlier_indx], iter_max_n=iter_max_n_full_loop, iter_convergence=iter_convergence)        
    
    return beta, outlier_indx, weights, residuals

def iwlls(measurements=None, design_matrix=None, weights=None, iter_max_n = 20, iter_convergence=1.0e-3):
    if measurements is None:
        return -1
    if design_matrix is None:
        return -1
    if weights is None:
        weights = np.ones(measurements.shape)

    w = np.ones(measurements.shape)*weights
    log_measurements = np.log(measurements)
    W = np.diag(w)
    # Initial Fit
    # (X'*W*X)\(X'*W)*log_measurements
    beta, residual, rank, sigma = linalg.lstsq( ((design_matrix.T).dot(W)).dot(design_matrix),
                                    ((design_matrix.T).dot(W)).dot(log_measurements), cond=None,
                                     overwrite_a=False, overwrite_b=False,
                                     check_finite=False, lapack_driver='gelsy')
    
    loop_criteria = True
    loop_n = 0
    while loop_criteria:
        loop_n += 1
        beta_previous = beta
        W = np.diag(w**2 * weights)
        beta, residual, rank, sigma = linalg.lstsq( ((design_matrix.T).dot(W)).dot(design_matrix),
                            ((design_matrix.T).dot(W)).dot(log_measurements), cond=None,
                                overwrite_a=False, overwrite_b=False,
                                check_finite=False, lapack_driver='gelsy')
        w = np.exp(design_matrix.dot(beta.T))
        loop_criteria = _check_loop_convergence(beta, beta_previous, loop_n, iter_convergence, iter_max_n)
    
    return beta

def _outlier_detection_calculate_residuals(measurements=None, design_matrix=None, 
            iter_max_n_inner_loop=5, iter_max_n_full_loop=20, iter_convergence=1.0e-3):
    # Iteratively re-Weighted Linear Least Squares estimator for DT/DK models
    # REKINDLE by Tax et al., MRM 2015

    w = np.ones(measurements.shape)
    log_measurements = np.log(measurements)
    W = np.diag(w)
    # (X'*W*X)\(X'*W)*log_measurements
    beta, residual, rank, sigma = linalg.lstsq( ((design_matrix.T).dot(W)).dot(design_matrix),
                                    ((design_matrix.T).dot(W)).dot(log_measurements), cond=None,
                                     overwrite_a=False, overwrite_b=False,
                                     check_finite=False, lapack_driver='gelsy')
    
    full_loop_criteria = True
    full_loop_n = 0
    while full_loop_criteria:
        full_loop_n += 1
        beta_previous = beta

        # Step 2: Compute a robust estimate for homoscedastic regression using IRLS.
        beta, w, loop_n = _iwlls_robust_homoscedastic_regression_estimate(beta=beta,
                                        design_matrix=design_matrix, log_measurements=log_measurements, 
                                        iter_max_n_inner_loop=iter_max_n_inner_loop, iter_convergence=iter_convergence)
        #  Step 3: Transform variables for heteroscedasticity
        estimated_log_measurements = design_matrix.dot(beta.T)
        heteroscedastic_log_measurements = log_measurements / np.exp(-estimated_log_measurements)
        design_matrix_star = design_matrix / np.exp(-estimated_log_measurements[:,None]) #[:,None] performs elementwise division for each column
        
        # Step 4: Initial LLS fit in * domain
        w = np.ones(measurements.shape)
        W = np.diag(w)
        beta, residual, rank, sigma = linalg.lstsq( ((design_matrix_star.T).dot(W)).dot(design_matrix_star),
                                        ((design_matrix_star.T).dot(W)).dot(heteroscedastic_log_measurements), cond=None,
                                        overwrite_a=False, overwrite_b=False,
                                        check_finite=False, lapack_driver='gelsy')
        # Step 5: Compute a robust estimate for homoscedastic regression using IRLS.
        beta, w, loop_n = _iwlls_robust_homoscedastic_regression_estimate(beta=beta,
                                        design_matrix=design_matrix_star, log_measurements=heteroscedastic_log_measurements,
                                        iter_max_n_inner_loop=iter_max_n_inner_loop, iter_convergence=iter_convergence)

        #  Step 6: Check convergence
        full_loop_criteria = _check_loop_convergence(beta, beta_previous, full_loop_n, iter_convergence, iter_max_n_full_loop)

    estimated_log_measurements = design_matrix_star.dot(beta.T)
    log_measurement_residuals = heteroscedastic_log_measurements - estimated_log_measurements

    return log_measurement_residuals

def _iwlls_robust_homoscedastic_regression_estimate(beta, design_matrix, log_measurements,
                                             iter_max_n_inner_loop, iter_convergence):
    loop_criteria = True
    loop_n = 0
    w = np.ones(log_measurements.shape)  
    while loop_criteria:        
        # a. Calculate the residuals e in the linear domain
        estimated_log_measurements = design_matrix.dot(beta.T)
        log_measurement_residuals = log_measurements - estimated_log_measurements
        # b. Obtain an estimate of the dispersion of the residuals with MAD.
        C = 1.4826*np.median(np.abs(log_measurement_residuals - np.median(log_measurement_residuals)))
        if C < iter_convergence:
            # This check avoids the problem of the unnecessary b0 downweighting.
            break

        loop_n += 1
        beta_previous = beta
        # c. Recompute the weights according to Eq. [13].
        w = (1 / (1 + (log_measurement_residuals/C)**2)**2)        
        W = np.diag(w)
        # d. Perform WLLS fit with new weights
        beta, residual, rank, sigma = linalg.lstsq( ((design_matrix.T).dot(W)).dot(design_matrix),
                                    ((design_matrix.T).dot(W)).dot(log_measurements), cond=None,
                                     overwrite_a=False, overwrite_b=False,
                                     check_finite=False, lapack_driver='gelsy')
        # e. Check convergence
        loop_criteria = _check_loop_convergence(beta, beta_previous, loop_n, iter_convergence, iter_max_n_inner_loop)

    return beta, w, loop_n

def _check_loop_convergence(beta, beta_previous, loop_n, iter_convergence, iter_max_n):
    if loop_n >= iter_max_n:
        return False
    if np.all(np.abs(beta - beta_previous) <= iter_convergence*np.maximum(np.abs(beta), np.abs(beta_previous))):
        return False
    return True
