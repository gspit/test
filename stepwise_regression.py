"""
Stepwise Linear Regression Analysis
==================================
Pure Python (NumPy + SciPy) implementation replicating SPSS output format.
Data: F:/article/Model DATA/Regression Tm.xlsx
Dependent Variable: Tm (Polymer Melting Temperature, Column C)
Independent Variables: Columns D-U (18 molecular descriptors)
Output: Stepwise_Regression_Report.htm, README_Stepwise_Regression.txt

Usage: python stepwise_regression.py

Requirements: numpy, pandas, scipy (all included in Anaconda/miniconda)
"""

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from scipy.linalg import lstsq
import zipfile
import xml.etree.ElementTree as ET
import os

def read_xlsx(path):
    """Read xlsx file without openpyxl - uses zipfile + xml."""
    with zipfile.ZipFile(path, 'r') as z:
        shared = []
        if 'xl/sharedStrings.xml' in z.namelist():
            ss_root = ET.fromstring(z.read('xl/sharedStrings.xml'))
            NS = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
            for si in ss_root.findall('.//ns:si', {'ns': NS}):
                t_nodes = si.findall('.//ns:t', {'ns': NS})
                shared.append(''.join(t.text or '' for t in t_nodes))
        sheet_root = ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
        NS = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
        rows = []
        for row_el in sheet_root.findall('.//ns:row', {'ns': NS}):
            row_data = []
            for cell in row_el.findall('ns:c', {'ns': NS}):
                t_attr = cell.get('t', '')
                v_el = cell.find('ns:v', {'ns': NS})
                if v_el is not None and v_el.text:
                    if t_attr == 's':
                        idx = int(v_el.text)
                        row_data.append(shared[idx] if idx < len(shared) else '')
                    else:
                        row_data.append(v_el.text)
                else:
                    row_data.append('')
            rows.append(row_data)
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows[1:], columns=rows[0])
    for col in df.columns:
        try:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        except:
            pass
    return df

def ols_reg(X, y):
    """Ordinary Least Squares regression. Returns dict with all statistics."""
    n_obs, k = X.shape
    Xc = np.column_stack([np.ones(n_obs), X])
    coeffs, *_ = lstsq(Xc, y, lapack_driver='gelsy')
    y_pred = Xc @ coeffs
    resid = y - y_pred
    rank = np.linalg.matrix_rank(Xc)
    df_mod = rank - 1
    df_res = n_obs - rank
    ss_tot = np.sum((y - y.mean())**2)
    ss_reg = np.sum((y_pred - y.mean())**2)
    ss_res = np.sum(resid**2)
    ms_reg = ss_reg / df_mod if df_mod > 0 else 0
    ms_res = ss_res / df_res if df_res > 0 else 0
    rsq = ss_reg / ss_tot if ss_tot > 0 else 0
    rsq_adj = 1 - (ss_res/df_res)/(ss_tot/(n_obs-1)) if df_res>0 and n_obs>1 else 0
    f_val = ms_reg/ms_res if ms_res>0 else 0
    f_p = 1-scipy_stats.f.cdf(f_val,df_mod,df_res) if ms_res>0 else 1
    try:
        cov = ms_res*np.linalg.inv(Xc.T@Xc)
        se = np.sqrt(np.diag(cov))
    except:
        se = np.full(k+1, np.nan)
    t_vals = coeffs/se
    t_pvals = np.array([2*(1-scipy_stats.t.cdf(abs(t),df_res)) for t in t_vals])
    return dict(params=coeffs, se=se, tvalues=t_vals, pvalues=t_pvals,
        rsquared=rsq, rsquared_adj=rsq_adj, fvalue=f_val, f_pvalue=f_p,
        mse_resid=ms_res, df_model=df_mod, df_resid=df_res,
        fittedvalues=y_pred, resid=resid, ss_reg=ss_reg, ss_res=ss_res)

def betas_calc(params, Xdf, y, names):
    """Compute standardized (Beta) coefficients."""
    y_std = y.std(ddof=1)
    return [params[i+1]*Xdf[names[i]].std(ddof=1)/y_std if y_std>0 else 0 for i in range(len(names))]

def vif_func(X, names):
    """Compute Variance Inflation Factor for each predictor."""
    result = []
    for i, name in enumerate(names):
        others = np.delete(X, i, axis=1)
        yv = X[:, i]
        try:
            coefs, *_ = lstsq(np.column_stack([np.ones(len(yv)), others]), yv, lapack_driver='gelsy')
            pred = np.column_stack([np.ones(len(yv)), others]) @ coefs
            ss_tot = np.sum((yv-yv.mean())**2)
            ss_res = np.sum((yv-pred)**2)
            r2 = 1-ss_res/ss_tot if ss_tot>0 else 0
            vif_val = 1/(1-r2) if r2<1 else 1000
        except:
            vif_val = 1.0
        result.append({'Variable': name, 'VIF': vif_val})
    return result

def partial_corr(x, resid):
    """Compute partial correlation between x and residuals."""
    mask = ~np.isnan(x) & ~np.isnan(resid)
    x_c = x[mask]-x[mask].mean()
    r_c = resid[mask]-resid[mask].mean()
    denom = np.sqrt(np.sum(x_c**2)*np.sum(r_c**2))
    return np.dot(x_c,r_c)/denom if denom>0 else 0.0

def calc_correlation_matrix(Xdf, y, names):
    """Calculate Pearson correlation matrix."""
    n_vars = len(names) + 1
    corr_matrix = np.zeros((n_vars, n_vars))
    all_data = np.column_stack([Xdf[names].values, y])
    for i in range(n_vars):
        for j in range(n_vars):
            if i == j:
                corr_matrix[i, j] = 1.0
            elif j < i:
                continue
            else:
                xi = all_data[:, i]
                xj = all_data[:, j]
                mask = ~np.isnan(xi) & ~np.isnan(xj)
                if mask.sum() > 2:
                    corr = np.corrcoef(xi[mask], xj[mask])[0, 1]
                    corr_matrix[i, j] = corr if not np.isnan(corr) else 0
                    corr_matrix[j, i] = corr_matrix[i, j]
    return corr_matrix

def durbin_watson(residuals):
    """Calculate Durbin-Watson statistic."""
    diff = np.diff(residuals)
    return np.sum(diff**2) / np.sum(residuals**2)

def stepwise(Xdf, y, preds, thr_in=0.05, thr_out=0.10):
    """
    Stepwise linear regression.
    Entry: F-test p-value < threshold_in (default 0.05)
    Removal: F-test p-value > threshold_out (default 0.10)
    """
    curr, rem, hist = [], list(preds), []
    step = 0
    while rem:
        step += 1; changed = False
        # Try to enter
        best_p, best_v = None, None
        for v in rem:
            Xt = Xdf[curr+[v]].values.astype(float)
            r = ols_reg(Xt, y)
            p_enter = r['pvalues'][-1]
            if p_enter < thr_in and (best_p is None or p_enter < best_p):
                best_p = p_enter; best_v = v
        if best_v:
            curr.append(best_v); rem.remove(best_v); changed = True
            hist.append(('IN', best_v, best_p, step))
        # Try to remove
        if curr:
            Xc = Xdf[curr].values.astype(float)
            r = ols_reg(Xc, y)
            pvs = r['pvalues'][1:]
            worst_p, worst_v = None, None
            for i, v in enumerate(curr):
                if pvs[i] > thr_out and (worst_p is None or pvs[i] > worst_p):
                    worst_p = pvs[i]; worst_v = v
            if worst_v:
                curr.remove(worst_v); rem.append(worst_v); changed = True
                hist.append(('OUT', worst_v, worst_p, step))
        if not changed: break
    return curr, hist

if __name__ == '__main__':
    print("Loading data...")
    DATA_PATH = 'F:/article/Regression Tm.xlsx'
    df_raw = read_xlsx(DATA_PATH)
    # Remove duplicate of first row at end
    if len(df_raw) > 1 and df_raw.iloc[0].equals(df_raw.iloc[-1]):
        df_raw = df_raw.iloc[:-1].reset_index(drop=True)
    # Drop rows with missing Tm or MolWt
    df = df_raw.dropna(subset=['Tm', 'MolWt']).reset_index(drop=True)
    print(f"N = {len(df)}")

    Y_COL = 'Tm'
    PREDICTOR_COLS = [
        'MolWt', 'LogP', 'TPSA', 'MolMR', 'NumRotatableBonds',
        'NumRigidBonds', 'NumBridgeheadAtoms', 'NumSpiroAtoms',
        'NumRings', 'NumAromaticRings', 'NumAliphaticRings',
        'RingDensity', 'NumHDonors', 'NumHAcceptors',
        'NumHeteroatoms', 'BalabanJ', 'Chi0n', 'Kappa3'
    ]
    available_preds = [c for c in PREDICTOR_COLS if c in df.columns]
    Y = df[Y_COL].values.astype(float)
    X_all_df = df[available_preds].copy()
    n = len(Y)

    print("Running stepwise regression...")
    selected_vars, step_history = stepwise(X_all_df, Y, available_preds)
    print(f"Selected {len(selected_vars)} variables")
    for action, var, pval, step in step_history:
        print(f"  Step {step}: {action:4s} {var} p={pval:.6f}")

    # Compute per-step model stats
    all_stats = []; curr_list = []
    for var in selected_vars:
        curr_list.append(var)
        Xs = X_all_df[curr_list].values.astype(float)
        r = ols_reg(Xs, Y)
        prev_r2 = all_stats[-1]['rsquared'] if all_stats else 0
        r2_chg = r['rsquared'] - prev_r2
        df2 = n - len(curr_list) - 1
        f_chg = (r2_chg/(1-r['rsquared']))*df2 if r2_chg>1e-10 and df2>0 else 0
        p_chg = 1-scipy_stats.f.cdf(f_chg,1,df2) if df2>0 else 1
        all_stats.append({'vars':list(curr_list),'R':np.sqrt(r['rsquared']),
            'rsquared':r['rsquared'],'rsquared_adj':r['rsquared_adj'],
            'SE':np.sqrt(r['mse_resid']),'F':r['fvalue'],'df_model':int(r['df_model']),
            'df_resid':int(r['df_resid']),'f_pvalue':r['f_pvalue'],
            'r2_chg':r2_chg,'f_chg':f_chg,'p_chg':p_chg,'r':r})

    final_r = all_stats[-1]['r']
    final_vars = selected_vars
    final_coeffs = final_r['params'][1:]
    intercept = final_r['params'][0]
    
    # Calculate additional statistics
    excluded_vars = [v for v in available_preds if v not in final_vars]
    
    # Descriptive statistics
    desc_stats = {}
    for col in [Y_COL] + final_vars:
        data = df[col].dropna()
        desc_stats[col] = {
            'mean': data.mean(),
            'std': data.std(ddof=1),
            'min': data.min(),
            'max': data.max(),
            'n': len(data)
        }
    
    # Correlation matrix for final variables
    corr_matrix = calc_correlation_matrix(X_all_df, Y, final_vars)
    
    # Residual analysis
    residuals = final_r['resid']
    y_pred = final_r['fittedvalues']
    standardized_resid = residuals / np.sqrt(final_r['mse_resid'])
    dw_stat = durbin_watson(residuals)
    
    # Partial correlations
    partial_corrs = {}
    for var in final_vars:
        idx = final_vars.index(var)
        X_other = np.delete(X_all_df[final_vars].values, idx, axis=1)
        r_other = ols_reg(X_other, Y)
        partial_corrs[var] = partial_corr(X_all_df[var].values, r_other['resid'])
    
    # VIF calculation
    vif_results = vif_func(X_all_df[final_vars].values, final_vars)
    
    # Save comprehensive report
    output_file = os.path.join(os.path.dirname(DATA_PATH), 'Stepwise_Regression_Full_Report.txt')
    with open(output_file, 'w', encoding='utf-8') as f:
        # Title
        f.write("="*90 + "\n")
        f.write("                    STEPWISE LINEAR REGRESSION ANALYSIS\n")
        f.write("="*90 + "\n\n")
        
        # Data Summary
        f.write("DATA SUMMARY\n")
        f.write("-"*90 + "\n")
        f.write(f"Data File: {DATA_PATH}\n")
        f.write(f"Dependent Variable: {Y_COL}\n")
        f.write(f"Total Observations: {n}\n")
        f.write(f"Independent Variables: {len(available_preds)}\n")
        f.write(f"Variables Entered: {len(final_vars)}\n")
        f.write(f"Variables Excluded: {len(excluded_vars)}\n\n")
        
        # Descriptive Statistics
        f.write("\nDESCRIPTIVE STATISTICS\n")
        f.write("-"*90 + "\n")
        f.write(f"{'Variable':<25} {'Mean':>12} {'Std. Deviation':>12} {'Minimum':>12} {'Maximum':>12} {'N':>8}\n")
        f.write("-"*90 + "\n")
        for col in [Y_COL] + final_vars:
            s = desc_stats[col]
            f.write(f"{col:<25} {s['mean']:>12.4f} {s['std']:>12.4f} {s['min']:>12.4f} {s['max']:>12.4f} {s['n']:>8}\n")
        f.write("-"*90 + "\n")
        
        # Correlation Matrix
        f.write("\nCORRELATION MATRIX\n")
        f.write("-"*90 + "\n")
        header = f"{'':<20}" + "".join([f"{v:>10}" for v in final_vars]) + f"{Y_COL:>10}\n"
        f.write(header)
        f.write("-"*90 + "\n")
        for i, var in enumerate(final_vars):
            line = f"{var:<20}"
            for j in range(len(final_vars)):
                line += f"{corr_matrix[i,j]:>10.3f}"
            line += f"{corr_matrix[len(final_vars),i]:>10.3f}\n"
            f.write(line)
        # Y row
        line = f"{Y_COL:<20}"
        for j in range(len(final_vars)):
            line += f"{corr_matrix[len(final_vars),j]:>10.3f}"
        line += f"{1.0:>10.3f}\n"
        f.write(line)
        f.write("-"*90 + "\n")
        
        # Variables Entered/Removed
        f.write("\nVARIABLES ENTERED/REMOVED\n")
        f.write("-"*90 + "\n")
        f.write(f"{'Step':>6} {'Action':<10} {'Variable':<25} {'p-value':>12}\n")
        f.write("-"*90 + "\n")
        for action, var, pval, step in step_history:
            f.write(f"{step:6d} {action:<10} {var:<25} {pval:>12.6f}\n")
        f.write("-"*90 + "\n")
        f.write(f"Method: Stepwise (criteria: probability-of-F-to-enter <= 0.050, probability-of-F-to-remove >= 0.100)\n\n")
        
        # Model Summary
        f.write("\nMODEL SUMMARY\n")
        f.write("-"*90 + "\n")
        f.write(f"{'Model':>6} {'R':>10} {'R Square':>10} {'Adj R Square':>14} {'Std. Error':>12} {'R Square Change':>16}\n")
        f.write("-"*90 + "\n")
        for i, stat in enumerate(all_stats, 1):
            f.write(f"{i:6d} {stat['R']:>10.4f} {stat['rsquared']:>10.4f} {stat['rsquared_adj']:>14.4f} ")
            f.write(f"{stat['SE']:>12.4f} {stat['r2_chg']:>16.4f}\n")
        f.write("-"*90 + "\n")
        
        # ANOVA
        f.write("\nANOVA\n")
        f.write("-"*90 + "\n")
        f.write(f"{'Model':>6} {'Sum of Squares':>16} {'df':>6} {'Mean Square':>14} {'F':>12} {'Sig.':>10}\n")
        f.write("-"*90 + "\n")
        for i, stat in enumerate(all_stats, 1):
            f.write(f"{i:6d} Regression  {stat['r']['ss_reg']:>14.4f} {stat['df_model']:>6d} {stat['r']['mse_resid']*stat['df_model']/stat['df_model']:>14.4f} ")
            f.write(f"{stat['F']:>12.4f} {stat['f_pvalue']:>10.6f}\n")
            f.write(f"      Residual    {stat['r']['ss_res']:>14.4f} {stat['df_resid']:>6d} {stat['r']['mse_resid']:>14.4f}\n")
            f.write(f"      Total       {stat['r']['ss_reg']+stat['r']['ss_res']:>14.4f} {n-1:>6d}\n")
        f.write("-"*90 + "\n")
        
        # Coefficients
        f.write("\nCOEFFICIENTS\n")
        f.write("-"*110 + "\n")
        f.write(f"{'Model':>6} {'':<20} {'Unstandardized':>16} {'Standardized':>14} {'t':>10} {'Sig.':>10} ")
        f.write(f"{'Collinearity':>14}\n")
        f.write(f"{'':<27} {'B':>10} {'Std. Error':>10} {'Beta':>10} {'':>10} {'':>10} {'VIF':>10}\n")
        f.write("-"*110 + "\n")
        
        betas = betas_calc(final_r['params'], X_all_df, Y, final_vars)
        for i, (var, coef, se, beta, t, p, vif) in enumerate(zip(final_vars, final_coeffs, 
                                                                 final_r['se'][1:], betas, 
                                                                 final_r['tvalues'][1:], 
                                                                 final_r['pvalues'][1:],
                                                                 [v['VIF'] for v in vif_results])):
            if i == 0:
                f.write(f"{len(all_stats):6d} (Constant)  {intercept:>10.4f} {final_r['se'][0]:>10.4f} {'':>10} {final_r['tvalues'][0]:>10.4f} {final_r['pvalues'][0]:>10.6f} {'':>14}\n")
            f.write(f"{'':6s} {var:<20} {coef:>10.4f} {se:>10.4f} {beta:>10.4f} {t:>10.4f} {p:>10.6f} {vif:>10.4f}\n")
        f.write("-"*110 + "\n")
        f.write(f"Dependent Variable: {Y_COL}\n\n")
        
        # Excluded Variables
        f.write("\nEXCLUDED VARIABLES\n")
        f.write("-"*110 + "\n")
        f.write(f"{'Model':>6} {'Variable':<25} {'Beta In':>10} {'t':>10} {'Sig.':>10} {'Partial Correlation':>20}\n")
        f.write("-"*110 + "\n")
        for var in excluded_vars:
            X_test = X_all_df[final_vars + [var]].values
            r_test = ols_reg(X_test, Y)
            beta_in = r_test['params'][-1] * X_all_df[var].std() / Y.std()
            t_in = r_test['tvalues'][-1]
            p_in = r_test['pvalues'][-1]
            partial = partial_corr(X_all_df[var].values, residuals)
            f.write(f"{len(all_stats):6d} {var:<25} {beta_in:>10.4f} {t_in:>10.4f} {p_in:>10.6f} {partial:>20.4f}\n")
        f.write("-"*110 + "\n")
        
        # Residual Statistics
        f.write("\nRESIDUAL STATISTICS\n")
        f.write("-"*90 + "\n")
        f.write(f"{'':<15} {'Minimum':>12} {'Maximum':>12} {'Mean':>12} {'Std. Deviation':>14} {'N':>8}\n")
        f.write("-"*90 + "\n")
        f.write(f"{'Predicted Value':<15} {y_pred.min():>12.4f} {y_pred.max():>12.4f} {y_pred.mean():>12.4f} {y_pred.std():>14.4f} {n:>8}\n")
        f.write(f"{'Residual':<15} {residuals.min():>12.4f} {residuals.max():>12.4f} {residuals.mean():>12.4f} {residuals.std():>14.4f} {n:>8}\n")
        f.write(f"{'Std. Predicted Value':<15} {(y_pred-y_pred.mean())/y_pred.std():>12.4f} {((y_pred-y_pred.mean())/y_pred.std()).max():>12.4f} {0:>12.4f} {1:>14.4f} {n:>8}\n")
        f.write(f"{'Std. Residual':<15} {standardized_resid.min():>12.4f} {standardized_resid.max():>12.4f} {standardized_resid.mean():>12.4f} {standardized_resid.std():>14.4f} {n:>8}\n")
        f.write("-"*90 + "\n")
        f.write(f"Durbin-Watson Statistic: {dw_stat:.4f}\n\n")
        
        # Additional Model Information
        f.write("\nADDITIONAL MODEL INFORMATION\n")
        f.write("-"*90 + "\n")
        f.write(f"Final Regression Equation:\n")
        f.write(f"{Y_COL} = {intercept:.4f}")
        for var, coef in zip(final_vars, final_coeffs):
            sign = '+' if coef >= 0 else '-'
            f.write(f"\n          {sign} {abs(coef):.4f} * {var}")
        f.write(f"\n\nModel Fit Statistics:\n")
        f.write(f"  R² = {final_r['rsquared']:.4f}\n")
        f.write(f"  Adjusted R² = {final_r['rsquared_adj']:.4f}\n")
        f.write(f"  F({final_r['df_model']}, {final_r['df_resid']}) = {final_r['fvalue']:.4f}, p < {final_r['f_pvalue']:.6f}\n")
        f.write(f"  RMSE = {np.sqrt(final_r['mse_resid']):.4f}\n")
        f.write(f"  AIC = {n*np.log(final_r['ss_res']/n) + 2*(len(final_vars)+1):.4f}\n")
        f.write(f"  BIC = {n*np.log(final_r['ss_res']/n) + np.log(n)*(len(final_vars)+1):.4f}\n\n")
        
        # Residual plots data (for reference)
        f.write("\nRESIDUAL PLOTS DATA\n")
        f.write("-"*90 + "\n")
        f.write("Histogram Data for Standardized Residuals (*zresid):\n")
        f.write(f"{'Bin Range':<20} {'Frequency':>10} {'Percent':>10}\n")
        f.write("-"*90 + "\n")
        hist_bins = [-3, -2, -1, 0, 1, 2, 3]
        freq, _ = np.histogram(standardized_resid, bins=hist_bins)
        for i in range(len(hist_bins)-1):
            pct = freq[i]/n*100
            f.write(f"{hist_bins[i]:>4.0f} to {hist_bins[i+1]:>4.0f}{'':<10} {freq[i]:>10d} {pct:>10.2f}\n")
        f.write("-"*90 + "\n\n")
        
        f.write("Note: For scatter plots, use the following variables:\n")
        f.write("  - ZRESID: Standardized residuals\n")
        f.write("  - ZPRED: Standardized predicted values\n")
        f.write("  - DRESID: Deleted residuals\n\n")
        
        f.write("="*90 + "\n")
        f.write("                         END OF REPORT\n")
        f.write("="*90 + "\n")
    
    print(f"\n✓ Full report saved to: {output_file}")
    print("\n=== Final Regression Equation ===")
    print(f"Tm = {intercept:.4f}")
    for var, coef in zip(final_vars, final_coeffs):
        sign = '+' if coef >= 0 else '-'
        print(f"    {sign} {abs(coef):.4f} * {var}")
    print(f"\nFinal: R2={final_r['rsquared']:.4f}, adjR2={final_r['rsquared_adj']:.4f}")
    print("Done! See output files.")
