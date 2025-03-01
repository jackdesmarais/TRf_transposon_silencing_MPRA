"""Library analysis module for MPRA data.

This module provides functions and classes for analyzing MPRA library data, including:
- Statistical analysis functions
- Plotting utilities 
- Library class for managing MPRA datasets
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib import rc_context
from matplotlib.colors import CenteredNorm

import statsmodels.api as sm
import seaborn as sns
from scipy.stats import pearsonr,gaussian_kde, barnard_exact

from math import ceil
from copy import copy

def count_over_thresh(x, threshold=0):
    """Count values over threshold in array.
    
    Args:
        x (array-like): Input array
        threshold (float, optional): Threshold value. Defaults to 0.
        
    Returns:
        array-like: Count of values over threshold. Returns NaN if all values are NaN.
    """
    over_thresh = x>threshold
    sums = over_thresh.sum()
    all_nans = x.isna().all()
    sums[all_nans] = np.nan
    return(sums)

def make_density(x):
    """Calculate kernel density estimate for array.
    
    Args:
        x (array-like): Input array
        
    Returns:
        pandas.Series: Density estimates for input values
    """
    nan_mask = np.logical_not(np.isnan(x))
    x = x[nan_mask]
    if len(x)<2:
        density = pd.Series(np.ones((len(x))))
    else:
        density = pd.Series(gaussian_kde(x).evaluate(x))
    density.index=x.index
    return(density)

def line_plotter(x, y, **kwargs):
    """Plot line through points.
    
    Args:
        x (array-like): x coordinates
        y (array-like): y coordinates
        **kwargs: Additional arguments passed to plt.axline()
    """
    nan_mask = np.logical_not(x.isna()|y.isna())
    if nan_mask.sum()>1:
        x = x[nan_mask]
        y=y[nan_mask]
        ax = plt.gca()
        ax.axline(**kwargs)

def correspondance_plotter(x,y,**kwargs):
    """Plot correlation between x and y values.
    
    Args:
        x (array-like): x coordinates
        y (array-like): y coordinates
        **kwargs: Additional arguments passed to plotting functions
    """
    nan_mask = np.logical_not(x.isna()|y.isna())
    if nan_mask.sum()>1:
        x = x[nan_mask]
        y=y[nan_mask]
        ax = plt.gca()

        r = pearsonr(x,y)
        plt.text(0.1,0.9,f'R:{r.statistic:.2f}',transform=ax.transAxes)

def density_scatter(x,y,**kwargs):
    """Create density scatter plot.
    
    Args:
        x (array-like): x coordinates
        y (array-like): y coordinates
        **kwargs: Additional arguments passed to sns.scatterplot()
    """
    nan_mask = np.logical_not(x.isna()|y.isna())
    if nan_mask.sum()>1:
        x = x[nan_mask]
        y=y[nan_mask]

        if nan_mask.sum()>2:
            d2 = np.stack([x,y])
            density = gaussian_kde(d2).evaluate(d2)
        else:
            density = np.ones(x.shape)

        sort_mask = np.argsort(density)

        sns.scatterplot(x=x.values[sort_mask], y=y.values[sort_mask], hue=density[sort_mask], **kwargs)

def OLS_plotter(x,y,ols_model=None, **kwargs):
    """Plot OLS regression line with confidence intervals.
    
    Args:
        x (array-like): x coordinates
        y (array-like): y coordinates
        ols_model (statsmodels.regression.linear_model.RegressionResultsWrapper, optional): Pre-fit OLS model
        **kwargs: Additional arguments passed to plt.plot()
    """
    nan_mask = np.logical_not(x.isna()|y.isna())
    if nan_mask.sum()>1:
        x = x[nan_mask]
        y=y[nan_mask]

        x_sorter= np.argsort(x.values)
        X = sm.add_constant(x.values[x_sorter])
        
       
        if ols_model is None:
            ols_model = sm.OLS(y.values[x_sorter], X).fit()

        pred = ols_model.get_prediction(X).summary_frame()
        plt.plot(x.values[x_sorter], pred['mean'], **kwargs)
        plt.plot(x.values[x_sorter], pred['mean_ci_lower'], linestyle='--', **kwargs)
        plt.plot(x.values[x_sorter], pred['mean_ci_upper'], linestyle='--', **kwargs)

def gini(sorted_arr):
    """Calculate Gini coefficient for sorted array.
    
    Args:
        sorted_arr (array-like): Sorted input array
        
    Returns:
        float: Gini coefficient
    """
    n = sorted_arr.size
    coef_ = 2. / n
    const_ = (n + 1.) / n
    weighted_sum = sum([(i+1)*yi for i, yi in enumerate(sorted_arr)])
    return coef_*weighted_sum/(sorted_arr.sum()) - const_

def multi_gini(xs,labels,title='Library skew', scatter_kwargs=dict()):
    """Plot Lorenz curves and calculate Gini coefficients for multiple arrays.
    
    Args:
        xs (list): List of arrays to analyze
        labels (list): Labels for each array
        title (str, optional): Plot title. Defaults to 'Library skew'.
        scatter_kwargs (dict, optional): Additional arguments passed to plt.scatter()
        
    Returns:
        tuple: (list of Gini coefficients, matplotlib figure)
    """
    kwargs = dict(marker='.')
    kwargs.update(scatter_kwargs)
    gs = []
    fig, ax = plt.subplots()
    ax.plot([0,1], [0,1], color='k',label='Perfectly uniform')
    for X_u, label in zip(xs, labels):
        X = X_u.copy()
        X = X[np.logical_not(np.isnan(X))]
        X.sort()
        X_lorenz = X.cumsum() / X.sum()
        X_lorenz = np.insert(X_lorenz, 0, 0)
        gs.append(gini(X))
        ax.scatter(np.arange(X_lorenz.size)/(X_lorenz.size-1), X_lorenz, 
            label=str(label)+' Gini: %.3f'%gs[-1],**kwargs)
    plt.xlabel('Sorted variant fraction')
    plt.ylabel('Cumulative read fraction')
    plt.legend(loc='upper left')
    plt.title(title)
    return(gs, fig)

def depth_plot(xs,labels,title='Library depth'):
    """Plot read depth distributions.
    
    Args:
        xs (list): List of arrays containing read counts
        labels (list): Labels for each array
        title (str, optional): Plot title. Defaults to 'Library depth'.
        
    Returns:
        tuple: (figure with read counts, figure with read fractions)
    """
    fig_read, ax_reads = plt.subplots()
    fig_frac, ax_fraction = plt.subplots()
    for X_u, label in zip(xs, labels):
        sort_obs = sorted(X_u[np.logical_not(np.isnan(X_u))],reverse=True)

        plt.sca(ax_reads)
        plt.loglog(np.arange(len(sort_obs))+1, sort_obs, '.', label=label)

        plt.sca(ax_fraction)
        plt.loglog(np.arange(len(sort_obs))+1, sort_obs/sum(sort_obs), '.', label=label)

    plt.sca(ax_reads)
    plt.xlabel('Variant rank-order')
    plt.ylabel('Variant read count')
    plt.legend(loc='lower left')
    plt.title(title)

    plt.sca(ax_fraction)
    plt.xlabel('Variant rank-order')
    plt.ylabel('Variant observation fraction')
    plt.legend(loc='lower left')
    plt.title(title)
    return(fig_read, fig_frac)

def make_ROC(true_class, metrics, labels, make_plots, plot_thresh=False):
    """Generate ROC curves and calculate AUC scores.
    
    Args:
        true_class (array-like): Binary class labels
        metrics (list): List of metric arrays to evaluate
        labels (list): Labels for each metric
        make_plots (bool): Whether to generate plots
        plot_thresh (bool, optional): Whether to plot threshold curves. Defaults to False.
        
    Returns:
        dict: AUC scores if make_plots=False
        tuple: (AUC scores, ROC figure, threshold figure, PRC figure, PRC threshold figure, histogram figures) if make_plots=True
    """
    if make_plots:
        print('Number of true positives: %s'%sum(true_class))
        fig_roc, ax1 = plt.subplots()
        fig_roc_thresh, ax2 = plt.subplots()
        fig3_PRC, ax3 = plt.subplots()
        fig4_PRC_Thresh, ax4 = plt.subplots()
        hist_figs, hist_axes = zip(*[plt.subplots() for i in range(len(metrics))])
    aucs = {}
    for i, (metric, label) in enumerate(zip(metrics,labels)):
        nan_mask = np.logical_not(np.isnan(metric))
        metric = metric[nan_mask]
        metric_classes = true_class[nan_mask]
        if make_plots:
            hist_ax = hist_axes[i]
            plt.sca(ax1)
        mask = np.argsort(metric)
        decision_curve = metric_classes[mask]
        tpr = np.cumsum(np.flip(decision_curve))/sum(metric_classes)
        fpr = np.cumsum(np.flip(decision_curve == 0))/sum(metric_classes == 0)
        auc = np.trapz(tpr,fpr)
        aucs[label+' ROC-AUC'] = auc
        if make_plots:
            if plot_thresh:
                plt.plot(fpr, tpr, label='%s Sorting AUC: %.2f'%(label, auc))
            else:
                plt.plot(fpr, tpr, label='%s AUC: %.2f'%(label, auc))

            thresholds = np.arange(np.nanmin(metric),np.nanmax(metric),(np.nanmax(metric)-np.nanmin(metric))/1000)[::-1]
            thresholds[0]=max(metric)

            tpr = np.array([sum(metric_classes[metric>=t])/sum(metric_classes) for t in thresholds])
            fpr = np.array([sum(metric_classes[metric>=t]== 0)/sum(metric_classes == 0) for t in thresholds])
            if plot_thresh:
                plt.plot(fpr, tpr, label='%s Thresh AUC: %.2f'%(label, np.trapz(tpr,fpr)))

            plt.sca(ax2)
            plt.plot(thresholds, tpr, label='%s True positive rate'%label)
            plt.plot(thresholds, fpr, label='%s False positive rate'%label)

            roc = np.subtract(tpr,fpr)
            t = thresholds[roc == np.nanmax(roc)]

            tmin = min(t)
            tmax = max(t)

            print('%s The optimal threshold range is from %s to %s'%(label, tmin, tmax))
            print('The ROCs at specific thresholds are- opt min: %s, opt max: %s'%(roc[thresholds==tmin],
                                                                                   roc[thresholds==tmax]))

            print('The TPRs at specific thresholds are- opt min: %s, opt max: %s'%(tpr[thresholds==tmin],
                                                                                   tpr[thresholds==tmax]))

            print('The FPRs at specific thresholds are- opt min: %s, opt max: %s'%(fpr[thresholds==tmin],
                                                                                   fpr[thresholds==tmax]))
            plt.sca(hist_ax)
            edges = np.histogram_bin_edges(metric,bins="auto")
            sns.histplot(metric, label = 'full distribution', 
                         bins=edges, kde=True)
            sns.histplot(metric[metric_classes], color='green', 
                         label = 'Sense mutants', bins=edges, kde=True)
            sns.histplot(metric[np.logical_not(metric_classes)], 
                         color='orange', label = 'Nonsense mutants', bins=edges, kde=True)

            plt.sca(ax2) 

            plt.sca(ax3)
        decision_curve = metric_classes[mask]
        tpr = np.cumsum(np.flip(decision_curve))/sum(metric_classes)
        precision = np.divide(np.cumsum(np.flip(decision_curve)), [i+1 for i in range(len(decision_curve))])
        auc = np.trapz(precision, tpr)
        aucs[label+' PRC-AUC'] = auc
        
        if make_plots:
            plt.plot(tpr, precision, label='%s Sorting AUC: %.2f'%(label, auc))

            tpr = np.array([sum(metric_classes[metric>=t])/sum(metric_classes) for t in thresholds])
            precision = np.array([sum(metric_classes[metric>=t])/sum(metric>=t) for t in thresholds])
            plt.plot(tpr, precision, label='%s Threshold AUC: %.2f'%(label, np.trapz(precision, tpr)))

            plt.sca(ax4)
            plt.plot(thresholds, tpr, label='%s True positive rate'%label)
            plt.plot(thresholds, precision, label='%s Precision'%label)

            f_score = np.divide(np.multiply(np.multiply(tpr,precision),2),np.add(tpr,precision))
            ft = thresholds[f_score == np.nanmax(f_score)]
            
            gm = np.sqrt(np.multiply(tpr,precision))
            gt = thresholds[gm == np.nanmax(gm)]
            
            
            tmin = min(t)
            tmax = max(t)

            
            plt.sca(hist_ax)
            tmin = min(ft)
            tmax = max(ft)
            print('The optimal f-score threshold range is from %s to %s'%(tmin, tmax))
            
            tmin = min(gt)
            tmax = max(gt)
            print('The optimal g-mean threshold range is from %s to %s'%(tmin, tmax))
                
            plt.title(label)
            plt.legend()

            plt.sca(ax4)
            if tmin == tmax:
                plt.axvline(tmin, c= 'r', alpha=0.5, label='optimal PRC thresh')
            elif tmin < tmax:
                plt.axvspan(tmin, tmax, facecolor='r', alpha=0.5, label='optimal PRC thresh')
            else:
                print('roc thresh error')
    
    if make_plots:
        print()
        print('Aggregate stats')
        plt.sca(ax1)
        plt.xlabel('False positive rate',fontsize=22)
        plt.ylabel('True positive rate',fontsize=22)
        plt.title('ROC',fontsize=22)
        plt.legend(loc='lower right', title='Area under the curve')

        plt.sca(ax2)
        plt.xlabel('Threshold',fontsize=22)
        plt.ylabel('Rate',fontsize=22)
        plt.title('ROC',fontsize=22)
        plt.legend(loc='upper right')

        plt.sca(ax3)
        plt.ylabel('Precision',fontsize=22)
        plt.xlabel('Recall',fontsize=22)
        plt.title('PRC',fontsize=22)
        plt.legend(loc='lower left', title='Area under the curve')

        plt.sca(ax4)
        plt.xlabel('Threshold',fontsize=22)
        plt.ylabel('Rate',fontsize=22)
        plt.title('PRC',fontsize=22)
        plt.legend(loc='lower left')
        
        return(aucs, fig_roc, fig_roc_thresh, fig3_PRC, fig4_PRC_Thresh, hist_figs)

    return(aucs)
    

def thresh_collapse(x, thresh=-0.5):
    """Calculate fraction of values above threshold.
    
    Args:
        x (array-like): Input array
        thresh (float, optional): Threshold value. Defaults to -0.5.
        
    Returns:
        float: Fraction of values above threshold
    """
    return((x>thresh).sum()/len(x))

class Library:
    """Class for analyzing MPRA library data.
    
    Attributes:
        data_df (pandas.DataFrame): Library data
        replicates (dict): Replicate information
        id_cols (list): Column names for variant IDs
        group_cols (list): Column names for grouping variants
        rate_method (str): Method for calculating rates
        alphabet (list): Valid characters in sequences
        WT_seq (str): Wild-type sequence
        fitness_measure (str): Column name for fitness values
    """

    def __init__(self, data_df, replicates, id_cols, group_cols, rate_method, alphabet, WT_seq, fitness_measure):
        """Initialize Library object.
        
        Args:
            data_df (pandas.DataFrame): Library data
            replicates (dict): Replicate information
            id_cols (list): Column names for variant IDs
            group_cols (list): Column names for grouping variants
            rate_method (str): Method for calculating rates
            alphabet (list): Valid characters in sequences
            WT_seq (str): Wild-type sequence
            fitness_measure (str): Column name for fitness values
        """
        self.data_df = data_df
        self.replicates = replicates
        self.id_cols = id_cols
        self.group_cols = group_cols
        self.rate_method = rate_method
        self.alphabet = alphabet
        self.WT_seq = WT_seq
        self.fitness_measure = fitness_measure

    @classmethod
    def build_from_reads(cls, data_file, replicates, id_cols, group_cols, rate_method, alphabet, WT_seq, process_call, name_regex, sheet_name=None):
        """Build Library object from read count data file.
        
        Args:
            data_file (str): Path to data file (.xlsx, .csv, or .tsv)
            replicates (dict): Replicate information
            id_cols (list): Column names for variant IDs
            group_cols (list): Column names for grouping variants
            rate_method (str): Method for calculating rates
            alphabet (list): Valid characters in sequences
            WT_seq (str): Wild-type sequence
            process_call (str): Column name for processing
            name_regex (str): Regex pattern for extracting names
            sheet_name (str, optional): Excel sheet name. Defaults to None.
            
        Returns:
            Library: New Library object
            
        Raises:
            ValueError: If data_file is not .xlsx, .csv or .tsv
        """
        if data_file.endswith('.xlsx'):
            data = pd.read_excel(data_file,sheet_name=sheet_name).dropna(axis=0,how='any')
        elif data_file.endswith('.csv'):
            data = pd.read_csv(data_file).dropna(axis=0,how='any')
        elif data_file.endswith('.tsv'):
            data = pd.read_csv(data_file,sep='\t').dropna(axis=0,how='any')
        else:
            raise ValueError('Data file must be an .xlsx, .tsv or .csv file')
        
        normalized_df = pd.DataFrame()
        for id_col in id_cols:
            normalized_df[('meta',id_col)] = data[id_col]

        for replicate, tps in replicates.items():
            for i, tp in enumerate(tps):
                normalized_df[(f'{replicate}_counts',tp)] = data[tp]
                normalized_df[(f'{replicate}_abundance',tp)] = data[tp]/data[tp].sum()
        normalized_df.columns = pd.MultiIndex.from_tuples(normalized_df.columns)

        for replicate, tps in replicates.items():
            if rate_method == 'l2fc':
                normalized_df[('FC',replicate)] = normalized_df[(f'{replicate}_abundance',tps[-1])]/normalized_df[(f'{replicate}_abundance',tps[0])]
                normalized_df[('l2fc',replicate)] = np.log2(normalized_df[('FC',replicate)])
            elif rate_method == 'l10fc':
                normalized_df[('FC',replicate)] = normalized_df[(f'{replicate}_abundance',tps[-1])]/normalized_df[(f'{replicate}_abundance',tps[0])]
                normalized_df[('l10fc',replicate)] = np.log10(normalized_df[('FC',replicate)])
        normalized_df = normalized_df.sort_index(axis=1)
        meta = normalized_df[('meta',process_call)].str.extract(name_regex)
        for col in meta.columns:
            normalized_df[('meta',col)] = meta[col]
        normalized_df = normalized_df.sort_index(axis=1)

        self = cls(normalized_df, replicates, id_cols, group_cols, rate_method, alphabet, WT_seq, rate_method)
        return(self)
    
    def save(self, filename):
        """Save Library object to a pickle file.
        
        Args:
            filename (str): Path to save the library object
            
        Raises:
            ValueError: If filename does not end in .pkl
        """
        if not filename.endswith('.pkl'):
            raise ValueError('Filename must end in .pkl')
            
        import pickle
        with open(filename, 'wb') as f:
            pickle.dump(self, f)
    
    @classmethod
    def load(cls, filename):
        """Load Library object from a pickle file.
        
        Args:
            filename (str): Path to the library pickle file
            
        Returns:
            Library: Loaded Library object
            
        Raises:
            ValueError: If filename does not end in .pkl
            ValueError: If loaded object is not a Library instance
        """
        if not filename.endswith('.pkl'):
            raise ValueError('Filename must end in .pkl')
            
        import pickle
        with open(filename, 'rb') as f:
            loaded = pickle.load(f)
            
        if not isinstance(loaded, cls):
            raise ValueError(f'Loaded object is not a {cls.__name__} instance')
            
        return loaded
    
    
    @classmethod
    def merge_by_overlap(cls, library_dest, library_source, dest_position_shift, source_position_shift, new_WT, position_col='Position', 
                         automatically_rescale_linearly=False, automatically_rescale_by_controls=False, 
                         control_col = 'AA', control_val = ['WT','*']):
        """Merge two Library objects by overlapping positions.
        
        Args:
            library_dest (Library): Destination library
            library_source (Library): Source library
            dest_position_shift (int): Position shift for destination library
            source_position_shift (int): Position shift for source library
            new_WT (str): New wild-type sequence
            position_col (str, optional): Column name for positions. Defaults to 'Position'.
            automatically_rescale_linearly (bool, optional): Whether to rescale linearly. Defaults to False.
            automatically_rescale_by_controls (bool, optional): Whether to rescale by controls. Defaults to False.
            control_col (str, optional): Column name for controls. Defaults to 'AA'.
            control_val (list, optional): Control values. Defaults to ['WT','*'].
            
        Returns:
            tuple: (Merged Library object, pre-rescaling figure, post-rescaling figure)
            
        Raises:
            AssertionError: If libraries are incompatible
        """
        assert len(library_dest.replicates.keys()&library_source.replicates.keys())==0, 'There are overlapping replicate names'
        merged_replicates = library_dest.replicates|library_source.replicates
        assert library_dest.id_cols == library_source.id_cols, 'Ensure ID columns are the same between libraries in both the library.id_cols and in the library.data_df'
        merged_id_cols = library_dest.id_cols
        assert library_dest.group_cols == library_source.group_cols, 'Ensure Grouping columns are the same between libraries in both the library.group_cols and in the library.data_df'
        merged_group_cols = library_dest.group_cols
        assert library_dest.rate_method == library_source.rate_method, 'Ensure rate methods are the same between libraries in both the library.rate_method and in the library.data_df'
        merged_rate_method = library_dest.rate_method
        assert library_dest.fitness_measure == library_source.fitness_measure, 'Ensure fitness measures are the same between libraries in both the library.fitness_measure and in the library.data_df'
        merged_fitness_measure = library_dest.fitness_measure
        assert (library_dest.alphabet == library_source.alphabet).all(), 'Ensure alphabets are the same between libraries in both the library.alphabet'
        merged_alphabet = library_dest.alphabet

        dest_df = library_dest.data_df.copy()
        wt_mask = dest_df[('meta',position_col)]==''
        dest_df[('meta',position_col)][~wt_mask] = dest_df[('meta',position_col)][~wt_mask].astype(int)+dest_position_shift

        source_df = library_source.data_df.copy()
        wt_mask = source_df[('meta',position_col)]==''
        source_df[('meta',position_col)][~wt_mask] = source_df[('meta',position_col)][~wt_mask].astype(int)+source_position_shift

        merged_df = pd.concat([dest_df,source_df],ignore_index=True).sort_index(axis=1)
        lib_merged = cls(merged_df, merged_replicates, merged_id_cols, merged_group_cols, merged_rate_method, merged_alphabet, new_WT, merged_fitness_measure)
        if automatically_rescale_linearly:
            lib_merged, fig_pre, fig_post = lib_merged.rescale_multiple_replicates_by_group(library_dest.replicates.keys(), library_source.replicates.keys())
        elif automatically_rescale_by_controls:
            lib_merged, fig_pre, fig_post = lib_merged.rescale_multiple_replicates_by_controls(library_dest.replicates.keys(), library_source.replicates.keys(),
                                                                            control_col, control_val)
        return(lib_merged, fig_pre, fig_post)
    
    @property
    def fitness_frame(self):
        """Get fitness data frame.
        
        Returns:
            pandas.DataFrame: Fitness values for each variant
        """
        df = self.data_df.set_index([('meta',col) for col in self.group_cols+self.id_cols])[self.fitness_measure].rename(mapper=lambda x: f'{self.fitness_measure} {x}', 
                axis='columns',
                level=0,
                inplace=False)
        df.index.names = self.group_cols+self.id_cols
        return(df.reset_index())
    
    @property
    def rep_average_frame(self):
        """Get replicate-averaged fitness data frame.
        
        Returns:
            pandas.DataFrame: Mean and std of fitness values across replicates
        """
        df = self.fitness_frame
        df = df.set_index(self.group_cols+self.id_cols).agg(['mean','std'],axis=1).rename(mapper=lambda x: f'{self.fitness_measure} {x}', 
                axis='columns',
                level=0,
                inplace=False)
        return(df.reset_index())
    
    @property
    def group_average_frame(self):
        """Get group-averaged fitness data frame.
        
        Returns:
            pandas.DataFrame: Mean and std of fitness values across groups
        """
        df = self.fitness_frame
        df = df.groupby(self.group_cols).agg(['mean','std'],axis=1).rename(mapper=lambda x: f'{self.fitness_measure} {x}', 
                axis='columns',
                level=1,
                inplace=False)
        return(df.reset_index())

    @property
    def total_average_frame(self):
        """Get total averaged fitness data frame.
        
        Returns:
            pandas.DataFrame: Mean fitness values across all replicates and groups
        """
        df = self.rep_average_frame
        df = df.groupby(self.group_cols)[f'{self.fitness_measure} mean'].mean().rename(f'{self.fitness_measure} mean')
        return(df.reset_index())
    
    def check_skew(self, sns_context='notebook', rc_params={}, title='Library skew'):
        """Check library skew using Lorenz curves.
        
        Args:
            sns_context (str, optional): Seaborn context. Defaults to 'notebook'.
            rc_params (dict, optional): RC parameters. Defaults to {}.
            title (str, optional): Plot title. Defaults to 'Library skew'.
            
        Returns:
            matplotlib.figure.Figure: Lorenz curve plot
        """
        with rc_context(rc_params), sns.plotting_context(sns_context):

            ginis, fig = multi_gini([self.data_df[(f'{key}_counts', tp)][~self.data_df[(f'{key}_counts', tp)].isna()].values for key, tps in self.replicates.items() for tp in tps],
                             [f'{key} {tp}' for key, tps in self.replicates.items() for tp in tps],
                             title=title)
            
            
        return(fig)
    
    def check_depth(self, sns_context='notebook', rc_params={}, title='Read depth'):
        """Check read depth distributions.
        
        Args:
            sns_context (str, optional): Seaborn context. Defaults to 'notebook'.
            rc_params (dict, optional): RC parameters. Defaults to {}.
            title (str, optional): Plot title. Defaults to 'Read depth'.
            
        Returns:
            tuple: (figure with read counts, figure with read fractions)
        """
        with rc_context(rc_params), sns.plotting_context(sns_context):

            fig_read, fig_frac = depth_plot([self.data_df[(f'{key}_counts', tp)][~self.data_df[(f'{key}_counts', tp)].isna()].values for key, tps in self.replicates.items() for tp in tps],
                             [f'{key} {tp}' for key, tps in self.replicates.items() for tp in tps],
                             title=title)
            
        return(fig_read, fig_frac)
    
    def check_barcode_distribution_skew(self,sns_context='notebook', threshold=0, rc_params={}, title='Read depth'):
        """Check skew in barcode distribution using Lorenz curves.
        
        Args:
            sns_context (str, optional): Seaborn context. Defaults to 'notebook'.
            threshold (int, optional): Minimum read count threshold. Defaults to 0.
            rc_params (dict, optional): RC parameters. Defaults to {}.
            title (str, optional): Plot title. Defaults to 'Read depth'.
            
        Returns:
            matplotlib.figure.Figure: Lorenz curve plot
        """
        with rc_context(rc_params), sns.plotting_context(sns_context):
            barcode_counts_df = self.data_df.groupby([('meta',col) for col in self.group_cols])[[(k+'_counts', col) for k, v in self.replicates.items() for col in v]].apply(count_over_thresh, threshold=threshold)
            ginis, fig = multi_gini([barcode_counts_df[(k+'_counts', col)][~barcode_counts_df[(f'{k}_counts', col)].isna()].values for k, v in self.replicates.items() for col in v],
                             [f'{key} {tp}' for key, tps in self.replicates.items() for tp in tps],
                             title=title)
            
        return(fig)

    def check_barcode_distribution_depth(self, sns_context='notebook', threshold=0, rc_params={}, title='Read depth'):
        """Check read depth distributions for barcodes.
        
        Args:
            sns_context (str, optional): Seaborn context. Defaults to 'notebook'.
            threshold (int, optional): Minimum read count threshold. Defaults to 0.
            rc_params (dict, optional): RC parameters. Defaults to {}.
            title (str, optional): Plot title. Defaults to 'Read depth'.
            
        Returns:
            tuple: (figure with read counts, figure with read fractions)
        """
        with rc_context(rc_params), sns.plotting_context(sns_context):
            barcode_counts_df = self.data_df.groupby([('meta',col) for col in self.group_cols])[[(k+'_counts', col) for k, v in self.replicates.items() for col in v]].apply(count_over_thresh, threshold=threshold)
            fig_read, fig_frac = depth_plot([barcode_counts_df[(k+'_counts', col)][~barcode_counts_df[(f'{k}_counts', col)].isna()].values for k, v in self.replicates.items() for col in v],
                             [f'{key} {tp}' for key, tps in self.replicates.items() for tp in tps],
                             title=title)
            
        return(fig_read, fig_frac)
    
    def check_controls(self, control_col, positive_vals, negative_vals, sns_context='notebook', rc_params={}, plt_type='strip'):
        """Check control distributions and generate ROC curves.
        
        Args:
            control_col (str): Column containing control labels
            positive_vals (list): Values indicating positive controls
            negative_vals (list): Values indicating negative controls 
            sns_context (str, optional): Seaborn context. Defaults to 'notebook'.
            rc_params (dict, optional): RC parameters. Defaults to {}.
            plt_type (str, optional): Plot type ('strip' or 'swarm'). Defaults to 'strip'.
            
        Returns:
            tuple: (replicates plot, averaged plot, ROC curve, ROC threshold plot, PRC curve, PRC threshold plot, histogram figures)
        """

        control_mask = self.data_df[('meta', control_col)].isin(positive_vals+negative_vals)
        positives = self.data_df[('meta', control_col)][control_mask].isin(positive_vals)

        to_plot = self.data_df[self.fitness_measure][control_mask]
        with rc_context(rc_params), sns.plotting_context(sns_context):
            
            
            to_plot = self.data_df[control_mask][self.fitness_measure]
            to_plot.columns.name = 'Replicate'
            to_plot.index.name = 'Variant number'
            to_plot = to_plot.stack().rename(self.fitness_measure).reset_index()
            
            control_labels = positives.rename('Control status').replace({True:'+',False:'-'})
            control_labels.index.name='Variant number'

            to_plot = to_plot.join(control_labels, on='Variant number')
                

            if plt_type == 'swarm':
                fig, replicates_cat_p = plt.subplots()
                replicates_cat_p = sns.swarmplot(data=to_plot, y=self.fitness_measure, x = 'Control status', size=1.5, order=['+','-'], hue='Replicate', dodge=True)
                l,h = replicates_cat_p.get_legend_handles_labels()
                sns.boxplot(data=to_plot, y=self.fitness_measure, x = 'Control status', order=['+','-'], boxprops=dict(alpha=0),showfliers=False,whis=0, hue='Replicate', dodge=True)
                plt.legend(l,h)
                
            elif plt_type == 'strip':
                fig, replicates_cat_p = plt.subplots()
                replicates_cat_p = sns.stripplot(data=to_plot, y=self.fitness_measure, x = 'Control status', size=1.5, order=['+','-'], hue='Replicate', dodge=True)
                l,h = replicates_cat_p.get_legend_handles_labels()
                sns.boxplot(data=to_plot, y=self.fitness_measure, x = 'Control status', order=['+','-'], boxprops=dict(alpha=0),showfliers=False,whis=0, hue='Replicate', dodge=True)
                plt.legend(l,h)
            
            to_plot = self.data_df[control_mask][self.fitness_measure].mean(axis=1).rename(self.fitness_measure+' replicate averaged')#
            to_plot.index.name = 'Variant number'
            to_plot = to_plot.reset_index()
            to_plot = to_plot.join(control_labels, on='Variant number')
            if plt_type == 'swarm':
                fig, averaged_cat_p = plt.subplots()
                averaged_cat_p = sns.swarmplot(data=to_plot, y=self.fitness_measure+' replicate averaged', x = 'Control status', size=1.5, order=['+','-'], hue='Control status')
                sns.boxplot(data=to_plot, y=self.fitness_measure+' replicate averaged', x = 'Control status', order=['+','-'], boxprops=dict(alpha=0),showfliers=False,whis=0, hue='Control status',dodge=False)
                
            if plt_type == 'strip':
                fig, averaged_cat_p = plt.subplots()
                to_plot['density'] = to_plot.groupby('Control status', group_keys=False)[self.fitness_measure+' replicate averaged'].apply(make_density)
                to_plot = to_plot.sort_values('density')
                averaged_cat_p = sns.stripplot(data=to_plot, y=self.fitness_measure+' replicate averaged', x = 'Control status', order=['+','-'], hue='density',legend=False)
                sns.boxplot(data=to_plot, y=self.fitness_measure+' replicate averaged', x = 'Control status', order=['+','-'], boxprops=dict(alpha=0),showfliers=False,whis=0)

            to_plot = self.data_df[control_mask][self.fitness_measure]
            to_plot.columns.name = 'Replicate'
            to_plot.index.name = 'Variant number'
            # to_plot = to_plot.stack().rename(self.fitness_measure).reset_index()
            
            control_labels = positives.rename('Control status').replace({True:'+',False:'-'})
            control_labels.index.name='Variant number'

            aucs, fig_roc, fig_roc_thresh, fig3_PRC, fig4_PRC_Thresh, hist_figs = make_ROC((control_labels=='+').values, 
                                                    [to_plot[col].values for col in to_plot.columns]+[to_plot.mean(axis=1).values], 
                                                    [col for col in to_plot.columns]+[self.fitness_measure+' averaged across reps'], 
                                                    True)
    
        return(replicates_cat_p, averaged_cat_p, fig_roc, fig_roc_thresh, fig3_PRC, fig4_PRC_Thresh, hist_figs)
    
    def make_replicate_checks(self, sns_context='notebook', rc_params={}, grid_kwargs=dict()):
        """Generate replicate-level fitness plots.
        
        Args:
            sns_context (str, optional): Seaborn context. Defaults to 'notebook'.
            rc_params (dict, optional): RC parameters. Defaults to {}.
            grid_kwargs (dict, optional): Additional arguments for the grid. Defaults to {}.
            
        Returns:
            tuple: (unaveraged plot, averaged plot)
        """
        with rc_context(rc_params), sns.plotting_context(sns_context):
            if len(self.replicates)==2:
                to_plot = self.data_df[self.fitness_measure].dropna(how='any',axis='index')
                xname = list(self.replicates.keys())[0]
                yname = list(self.replicates.keys())[1]

                unaveraged_g = sns.JointGrid(data = to_plot, x=xname, y=yname, **grid_kwargs) 
                unaveraged_g.plot_joint(density_scatter,edgecolor="none", legend=False)
                unaveraged_g.plot_joint(line_plotter, xy1=(0,0), slope=1, linestyle='--',color='k', alpha=0.3)
                unaveraged_g.plot_joint(correspondance_plotter)
                unaveraged_g.plot_marginals(sns.histplot)   
                unaveraged_g.plot_joint(OLS_plotter, color='C0', alpha=0.3)
            elif len(self.replicates)>2:
                to_plot = self.data_df[self.fitness_measure]
                unaveraged_g = sns.PairGrid(to_plot, corner=True, **grid_kwargs)
                unaveraged_g.map_diag(sns.histplot)
                unaveraged_g.map_lower(density_scatter,edgecolor="none", legend=False)
                unaveraged_g.map_lower(correspondance_plotter)
                unaveraged_g.map_lower(line_plotter, xy1=(0,0), slope=1, linestyle='--',color='k', alpha=0.3)
                unaveraged_g.map_lower(OLS_plotter, color='C0', alpha=0.3)
                
            if self.group_cols is not None:
                if len(self.replicates)==2:
                    to_plot = self.data_df.set_index([('meta',col) for col in self.group_cols])[self.fitness_measure].reset_index().groupby([('meta',col) for col in self.group_cols]).mean()
                    to_plot.index.names = self.group_cols

                    to_plot = to_plot.add_suffix(' group averaged')

                    xname = list(self.replicates.keys())[0]+' group averaged'
                    yname = list(self.replicates.keys())[1]+' group averaged'
                    
                    averaged_g = sns.JointGrid(data = to_plot, x=xname, y=yname, **grid_kwargs) 
                    averaged_g.plot_joint(density_scatter,edgecolor="none", legend=False)
                    averaged_g.plot_joint(line_plotter, xy1=(0,0), slope=1, linestyle='--',color='k', alpha=0.3)
                    averaged_g.plot_joint(correspondance_plotter)
                    averaged_g.plot_marginals(sns.histplot)   
                    averaged_g.plot_joint(OLS_plotter, color='C0', alpha=0.3)
                    
                elif len(self.replicates)>2:
                    to_plot = self.data_df.set_index([('meta',col) for col in self.group_cols])[self.fitness_measure].reset_index().groupby([('meta',col) for col in self.group_cols]).mean()
                    to_plot.index.names = self.group_cols
                    to_plot = to_plot.add_suffix(' group averaged')

                    averaged_g = sns.PairGrid(to_plot, corner=True, **grid_kwargs)
                    averaged_g.map_diag(sns.histplot)
                    averaged_g.map_lower(density_scatter,edgecolor="none", legend=False)
                    averaged_g.map_lower(correspondance_plotter)
                    averaged_g.map_lower(line_plotter, xy1=(0,0), slope=1, linestyle='--',color='k', alpha=0.3)
                    averaged_g.map_lower(OLS_plotter, color='C0', alpha=0.3)
            
        return(unaveraged_g,averaged_g)
    
    def get_rescaler_by_replicate(self, dest_replicate, source_replicate):
        """Get a rescaling model for a specific replicate pair.
        
        Args:
            dest_replicate (str): Destination replicate name
            source_replicate (str): Source replicate name
            
        Returns:
            tuple: (OLS model, plot)
        """
        to_fit = self.data_df[self.fitness_measure]
        

        X = sm.add_constant(to_fit[source_replicate].values)
        ols_model = sm.OLS(to_fit[dest_replicate].values, X).fit()

        g = sns.JointGrid(data = to_fit, x=source_replicate, y=dest_replicate) 
        g.plot_joint(density_scatter,edgecolor="none", legend=False)
        g.plot_joint(line_plotter, xy1=(0,0), slope=1, linestyle='--',color='k', alpha=0.3)
        g.plot_joint(correspondance_plotter)
        g.plot_marginals(sns.histplot)   
        g.plot_joint(OLS_plotter, ols_model=ols_model, color='C0', alpha=0.3)
        return(ols_model,g)
    
    def rescale_by_replicate(self, dest_replicate, source_replicate):
        """Rescale data by a specific replicate pair.
        
        Args:
            dest_replicate (str): Destination replicate name
            source_replicate (str): Source replicate name
            
        Returns:
            tuple: (rescaled library, pre-rescaling plot, post-rescaling plot)
        """                
        ols_model, pre_fig = self.get_rescaler_by_replicate(dest_replicate, source_replicate)
        new_library = copy(self)
        rescaled_df = self.data_df.copy()
        rescaled_df[(new_library.fitness_measure,source_replicate)] = ols_model.predict(sm.add_constant(rescaled_df[[(new_library.fitness_measure,source_replicate)]].values))
        new_library.data_df = rescaled_df
        
        _, post_fig = new_library.get_rescaler_by_replicate(dest_replicate, source_replicate)
        return(new_library,pre_fig,post_fig)
    
    def get_rescaler_by_replicate_groups(self, dest_replicate, source_replicate):
        """Get a rescaling model for a specific replicate pair across groups.
        
        Args:
            dest_replicate (str): Destination replicate name
            source_replicate (str): Source replicate name
            
        Returns:
            tuple: (OLS model, plot)
        """
        to_fit = self.data_df.set_index([('meta',col) for col in self.group_cols])[self.fitness_measure]
        to_fit.index.names = self.group_cols
        to_fit = to_fit.reset_index().groupby(self.group_cols, group_keys=True).mean()
        
        X = sm.add_constant(to_fit[source_replicate].values)
        ols_model = sm.OLS(to_fit[dest_replicate].values, X).fit()
        

        g = sns.JointGrid(data = to_fit, x=source_replicate, y=dest_replicate) 
        g.plot_joint(density_scatter,edgecolor="none", legend=False)
        g.plot_joint(line_plotter, xy1=(0,0), slope=1, linestyle='--',color='k', alpha=0.3)
        g.plot_joint(correspondance_plotter)
        g.plot_marginals(sns.histplot)   
        g.plot_joint(OLS_plotter, ols_model=ols_model, color='C0', alpha=0.3)
        return(ols_model, g)
    
    def rescale_by_replicate_group(self, dest_replicate, source_replicate):
        """Rescale data by a specific replicate pair across groups.
        
        Args:
            dest_replicate (str): Destination replicate name
            source_replicate (str): Source replicate name
        
        Returns:
            tuple: (rescaled library, pre-rescaling plot, post-rescaling plot)
        """            
        ols_model, pre_fig = self.get_rescaler_by_replicate_groups(dest_replicate, source_replicate)
        new_library = copy(self)
        rescaled_df = self.data_df.copy()
        rescaled_df[(new_library.fitness_measure,source_replicate)] = ols_model.predict(sm.add_constant(rescaled_df[[(new_library.fitness_measure,source_replicate)]].values))
        new_library.data_df = rescaled_df
        
        _, post_fig = new_library.get_rescaler_by_replicate_groups(dest_replicate, source_replicate)
        return(new_library, pre_fig,post_fig)
    
    def get_multi_replicate_rescaler_by_groups(self, dest_replicates, source_replicates):
        """Get a rescaling model for multiple replicate pairs across groups.
        
        Args:
            dest_replicates (list): List of destination replicate names
            source_replicates (list): List of source replicate names
        
        Returns:
            tuple: (OLS model, plot)
        """
        to_rescale = self.data_df.set_index([('meta',col) for col in self.group_cols])[self.fitness_measure]
        dest_to_scale = to_rescale[dest_replicates]
        source_to_scale = to_rescale[source_replicates]

        dest_name = 'Dest replicates '+self.fitness_measure
        dest_to_fit = dest_to_scale.mean(axis=1).rename(dest_name)
        dest_to_fit.index.names = self.group_cols
        dest_to_fit = dest_to_fit.reset_index()
        dest_to_fit = dest_to_fit.groupby(self.group_cols, group_keys=True).mean()
        dest_to_fit = dest_to_fit.sort_index(axis=1)

        source_name = 'Source replicates '+self.fitness_measure
        source_to_fit = source_to_scale.mean(axis=1).rename(source_name)
        source_to_fit.index.names = self.group_cols
        source_to_fit = source_to_fit.reset_index()
        source_to_fit = source_to_fit.groupby(self.group_cols, group_keys=True).mean()
        source_to_fit = source_to_fit.sort_index(axis=1)

        to_fit = source_to_fit.join(dest_to_fit).dropna(how='any')
        

        X = sm.add_constant(to_fit[source_name].values)
        ols_model = sm.OLS(to_fit[dest_name].values, X).fit()

        g = sns.JointGrid(data = to_fit, x=source_name, y=dest_name) 
        g.plot_joint(density_scatter,edgecolor="none", legend=False)
        g.plot_joint(line_plotter, xy1=(0,0), slope=1, linestyle='--',color='k', alpha=0.3)
        g.plot_joint(correspondance_plotter)
        g.plot_marginals(sns.histplot)   
        g.plot_joint(OLS_plotter, ols_model=ols_model, color='C0', alpha=0.3)
        
        return(ols_model, g)
    
    def rescale_multiple_replicates_by_group(self, dest_replicates, source_replicates,):
        """Rescale data by multiple replicate pairs across groups.
        
        Args:
            dest_replicates (list): List of destination replicate names
            source_replicates (list): List of source replicate names
            
        Returns:
            tuple: A tuple containing:
                - new_library (Library): A new Library object with rescaled data
                - pre_fig (matplotlib.figure): Plot showing fit before rescaling
                - post_fig (matplotlib.figure): Plot showing fit after rescaling
        """
        ols_model, pre_fig = self.get_multi_replicate_rescaler_by_groups(dest_replicates, source_replicates)
        new_library = copy(self)
        rescaled_df = self.data_df.copy()
        for rep in source_replicates:
            rescaled_df[(new_library.fitness_measure,rep)] = ols_model.predict(sm.add_constant(rescaled_df[[(new_library.fitness_measure,rep)]].values))
        new_library.data_df = rescaled_df
        _, post_fig = new_library.get_multi_replicate_rescaler_by_groups(dest_replicates, source_replicates,)
        return(new_library,pre_fig,post_fig)
    
    def get_multi_replicate_rescaler_by_controls(self, dest_replicates, source_replicates,control_col, control_values):
        """Get a rescaling model for multiple replicate pairs using control values.
        
        Args:
            dest_replicates (list): List of destination replicate names
            source_replicates (list): List of source replicate names
            control_col (str): Name of column containing control values
            control_values (list): List of values to use as controls
            
        Returns:
            tuple: A tuple containing:
                - ols_model (statsmodels.regression.linear_model.RegressionResultsWrapper): Fitted OLS model
                - g (seaborn.JointGrid): Plot showing the fit
        """
        to_rescale = self.data_df[self.data_df[('meta',control_col)].isin(control_values)].set_index(('meta',control_col))[self.fitness_measure]
        to_rescale.index.name = control_col

        dest_name = self.fitness_measure+' Dest replicates'
        source_name = self.fitness_measure+' Source replicates'

        dest_to_scale = to_rescale[dest_replicates].mean(axis=1).rename(dest_name).reset_index().groupby('AA').mean()
        source_to_scale = to_rescale[source_replicates].mean(axis=1).rename(source_name).reset_index().groupby('AA').mean()

        to_rescale = source_to_scale.join(dest_to_scale)

        X = sm.add_constant(to_rescale[source_name].values)
        ols_model = sm.OLS(to_rescale[dest_name].values, X).fit()

        g = sns.JointGrid(data = to_rescale, x=source_name, y=dest_name) 
        g.plot_joint(density_scatter,edgecolor="none", legend=False)
        g.plot_joint(line_plotter, xy1=(0,0), slope=1, linestyle='--',color='k', alpha=0.3)
        g.plot_joint(correspondance_plotter)
        g.plot_marginals(sns.histplot)   
        g.plot_joint(OLS_plotter, ols_model=ols_model, color='C0', alpha=0.3)
        
        return(ols_model, g)
    
    def rescale_multiple_replicates_by_controls(self, dest_replicates, source_replicates, control_col, control_values):
        """Rescale data by multiple replicate pairs using control values.
        
        Args:
            dest_replicates (list): List of destination replicate names
            source_replicates (list): List of source replicate names
            control_col (str): Name of column containing control values
            control_values (list): List of values to use as controls
            
        Returns:
            tuple: A tuple containing:
                - new_library (Library): A new Library object with rescaled data
                - pre_fig (matplotlib.figure): Plot showing fit before rescaling
                - post_fig (matplotlib.figure): Plot showing fit after rescaling
        """
        ols_model, pre_fig = self.get_multi_replicate_rescaler_by_controls(dest_replicates, source_replicates,control_col, control_values)
        new_library = copy(self)
        rescaled_df = self.data_df.copy()
        for rep in source_replicates:
            rescaled_df[(new_library.fitness_measure,rep)] = ols_model.predict(sm.add_constant(rescaled_df[[(new_library.fitness_measure,rep)]].values))
        new_library.data_df = rescaled_df
        _, post_fig = new_library.get_multi_replicate_rescaler_by_controls(dest_replicates, source_replicates,control_col, control_values)
        return(new_library, pre_fig, post_fig)


    def make_initial_skew_checks(self, sns_context='notebook', rc_params={}):
        """Make plots checking initial abundance vs fitness measure skew.
        
        Args:
            sns_context (str, optional): Seaborn plotting context. Defaults to 'notebook'.
            rc_params (dict, optional): Dictionary of rc parameters to use. Defaults to {}.
            
        Returns:
            tuple: A tuple containing:
                - self (Library): The Library object
                - figs (list): List of seaborn.JointGrid plots showing skew
        """
        with rc_context(rc_params), sns.plotting_context(sns_context):
            figs = []
            for key, tps in self.replicates.items():
                to_plot = self.data_df[[(f'{key}_abundance', tps[0]),(self.fitness_measure, key)]].dropna(how='any',axis='index')
                to_plot.columns = [f'{tps[0]} abundance', f'{key} {self.fitness_measure}']

                g = sns.JointGrid(data = to_plot, x=f'{tps[0]} abundance', y=f'{key} {self.fitness_measure}') 
                g.plot_joint(density_scatter, edgecolor="none", legend=False)
                # g.plot_joint(line_plotter, xy1=(0,0), slope=1, linestyle='--',color='k', alpha=0.3)
                g.plot_joint(correspondance_plotter)
                g.plot_marginals(sns.histplot) 
                figs.append(g)

        return(self, figs)

    def get_single_mutant_mutability_score(self, control_col,WT_vals, stop_vals, position_col, thresh, invert=False):
        """Calculate mutability scores for single mutants.
        
        Args:
            control_col (str): Name of column containing control values
            WT_vals (list): List of values indicating wild type
            stop_vals (list): List of values indicating stop codons
            position_col (str): Name of column containing position values
            thresh (float): Threshold value for mutability calculation
            invert (bool, optional): Whether to invert the mutability score. Defaults to False.
            
        Returns:
            pandas.DataFrame: DataFrame containing mutability scores by position
        """
        to_collapse = self.total_average_frame
        neg = to_collapse[control_col].isin(stop_vals)
        pos = to_collapse[control_col].isin(WT_vals)
        to_collapse = to_collapse[(~neg)&(~pos)]
        to_collapse = to_collapse.groupby(position_col, group_keys=True)[self.fitness_measure+' mean'].apply(thresh_collapse, thresh=thresh).rename('Mutability')
        if invert:
            to_collapse = 1-to_collapse
        return(to_collapse.reset_index())
    
    def get_single_mutant_property_enrichments(self, property_sets, thresh, position_col='Position', char_col='AA'):
        """Calculate property enrichments for single mutants.
        
        Args:
            property_sets (dict): Dictionary mapping property names to lists of amino acids
            thresh (float): Threshold value for enrichment calculation
            position_col (str, optional): Name of column containing position values. Defaults to 'Position'.
            char_col (str, optional): Name of column containing amino acid values. Defaults to 'AA'.
            
        Returns:
            pandas.DataFrame: DataFrame containing enrichment statistics by position and property
        """
        frame = self.total_average_frame
        checks = {}
        for n,s in property_sets.items():
            to_stat = frame[frame.AA!='*']
            to_stat = pd.crosstab((to_stat[f'{self.fitness_measure} mean']>thresh).rename('OverThresh'),[to_stat[position_col], frame[char_col].isin(s).rename('Property')])
            to_stat = to_stat.stack(position_col).unstack('OverThresh').stack('OverThresh').drop(index=[''])
            to_stat = to_stat.groupby(position_col).apply(lambda x: barnard_exact(x, pooled=False)).apply(lambda x: pd.Series((x.statistic, x.pvalue))).rename(columns={0:'Statistic',1:'P-Val'})
            checks[n]=to_stat

        return(pd.concat(checks))
    
    def draw_singles_heatmap(self, ax, cax, x_ax, mut_num_col, control_col, WT_vals, stop_vals, position_col, grid, cmap, norm,
                             seq_start, tick_interval, plot_pos):
        """Draw heatmap of single mutant effects.
        
        Args:
            ax (matplotlib.axes.Axes): Axes to draw heatmap on
            cax (matplotlib.axes.Axes): Axes for colorbar
            x_ax (matplotlib.axes.Axes): Axes for x-axis labels
            mut_num_col (str): Name of column containing number of mutations
            control_col (str): Name of column containing control values
            WT_vals (list): List of values indicating wild type
            stop_vals (list): List of values indicating stop codons
            position_col (str): Name of column containing position values
            grid (bool): Whether to draw grid lines
            cmap (matplotlib.colors.Colormap): Colormap to use
            norm (matplotlib.colors.Normalize): Normalization to use
            seq_start (int): Starting position in sequence
            tick_interval (int): Interval between ticks
            plot_pos (tuple): (start, end) positions to plot
        """
        singles_df = self.data_df[(self.data_df[('meta',mut_num_col)]=='1')|(self.data_df[('meta',control_col)].isin(WT_vals+stop_vals))]

        rep_averaged_df = singles_df.set_index([('meta',col) for col in self.group_cols])[self.fitness_measure].mean(axis=1).rename(self.fitness_measure)
        rep_averaged_df.index.names = self.group_cols
        rep_averaged_df = rep_averaged_df.reset_index()
        total_averaged_df = rep_averaged_df.groupby(self.group_cols, group_keys=True).agg(['mean','std'])
        total_averaged_df = total_averaged_df.sort_index(axis=1)

        wts = total_averaged_df.index.get_level_values(control_col).isin(WT_vals)
        
        total_averaged_df[(self.fitness_measure,'mean')] = total_averaged_df-total_averaged_df[wts][(self.fitness_measure,'mean')].mean()
        total_averaged_df = total_averaged_df[~wts]

        
        aa_encoding = {v:k for k,v in enumerate(self.alphabet)}
        aa = np.array([aa_encoding[aa] for aa in total_averaged_df.reset_index().AA.to_numpy()])
        wt_oh = np.array(list(self.WT_seq))[None,:]==self.alphabet[:,None]
        wt_aa_loc, wt_pos_loc = np.where(wt_oh)

        fitnesses = np.zeros((len(self.alphabet),len(self.WT_seq)))
        fitnesses[:,:]=np.nan
        fitnesses[aa.astype(int), total_averaged_df.reset_index()[position_col].to_numpy().astype(int)-1] = total_averaged_df[(self.fitness_measure,'mean')]
        fitnesses[wt_aa_loc,wt_pos_loc] = 0

        plt.sca(ax)
        fitnesses = np.zeros((len(self.alphabet),len(self.WT_seq)))
        fitnesses[:,:]=np.nan
        fitnesses[aa.astype(int), total_averaged_df.reset_index()[position_col].to_numpy().astype(int)-1] = total_averaged_df[(self.fitness_measure,'mean')]
        fitnesses[wt_aa_loc,wt_pos_loc] = 0

        if norm is None:
            norm=CenteredNorm(vcenter=0)
        if grid:
            plt.pcolormesh(fitnesses,shading='flat', edgecolors='k',
                    cmap=cmap,norm=norm, linewidth=grid)
            
            scatter_offset=-0.5
            tick_offset=0.5
            ax.invert_yaxis()
            # axes[1,0].grid(True, c='k')
        else:
            plt.imshow(fitnesses,interpolation='none', aspect='equal',
                    cmap=cmap,norm=norm)
            scatter_offset=0
            tick_offset=0

        cbar = plt.colorbar(orientation='vertical',cax=cax)
        cax.set_yscale('linear')
        cbar.set_label('Log$_{2}$ fold Change relative to wild type')

        
        plt.scatter(wt_pos_loc-scatter_offset, wt_aa_loc-scatter_offset, s=2, c='k')
        plt.yticks(np.arange(len(self.alphabet))-scatter_offset,self.alphabet)
        plt.xticks(np.arange(len(self.WT_seq))-scatter_offset,list(self.WT_seq))

        plt.xlabel('WT seq')

        ax.set_xlim(plot_pos)


        x_ax.patch.set_visible(False)
        x_ax.yaxis.set_visible(False)
        x_ax.spines[['right', 'top', 'left']].set_visible(False)
        ticks = np.arange(tick_interval*((seq_start+tick_interval)//tick_interval),seq_start+len(self.WT_seq),tick_interval)
        x_ax.set_xticks(ticks-seq_start+tick_offset, ticks)
        x_ax.set_xlim(plot_pos)
        x_ax.set_xlabel('Position')
        
        
    def draw_mutability_score(self, ax, cax, control_col, WT_vals, stop_vals, position_col, thresh,invert_mutability, grid, cmap, norm,plot_pos):
        """Draw heatmap of mutability scores.
        
        Args:
            ax (matplotlib.axes.Axes): Axes to draw heatmap on
            cax (matplotlib.axes.Axes): Axes for colorbar
            control_col (str): Name of column containing control values
            WT_vals (list): List of values indicating wild type
            stop_vals (list): List of values indicating stop codons
            position_col (str): Name of column containing position values
            thresh (float): Threshold value for mutability calculation
            invert_mutability (bool): Whether to invert mutability scores
            grid (bool): Whether to draw grid lines
            cmap (matplotlib.colors.Colormap): Colormap to use
            norm (matplotlib.colors.Normalize): Normalization to use
            plot_pos (tuple): (start, end) positions to plot
        """
        plt.sca(ax)
        to_collapse = self.get_single_mutant_mutability_score(control_col, WT_vals, stop_vals, position_col, thresh, invert=invert_mutability)
        properties = np.zeros((1,len(self.WT_seq)))
        properties[:,:]=np.nan
        properties[np.zeros(len(to_collapse[position_col].values)).astype(int),to_collapse[position_col].values.astype(int)-1] = to_collapse.Mutability.values
        if grid:
            plt.pcolormesh(properties,shading='flat', edgecolors='k',
                    cmap=cmap, linewidth=grid, norm=norm)
            
            # axes[0,0].set_aspect('equal')
        else:
            plt.imshow(properties,interpolation='none', 
                    cmap=cmap, norm=norm)
        ax.yaxis.set_visible(False)
        ax.xaxis.set_visible(False)
        cbar = plt.colorbar(orientation='vertical',cax=cax)
        cax.set_yscale('linear')
        ax.set_xlim(plot_pos)

    def draw_property_scores(self, ax, cax, property_sets, thresh, position_col, control_col, cmap, grid, norm, plot_pos):
        """Draw heatmap of property enrichment scores.
        
        Args:
            ax (matplotlib.axes.Axes): Axes to draw heatmap on
            cax (matplotlib.axes.Axes): Axes for colorbar
            property_sets (dict): Dictionary mapping property names to lists of amino acids
            thresh (float): Threshold value for enrichment calculation
            position_col (str): Name of column containing position values
            control_col (str): Name of column containing control values
            cmap (matplotlib.colors.Colormap): Colormap to use
            grid (bool): Whether to draw grid lines
            norm (matplotlib.colors.Normalize): Normalization to use
            plot_pos (tuple): (start, end) positions to plot
        """
        plt.sca(ax)
        to_collapse = self.get_single_mutant_property_enrichments(property_sets, thresh, position_col=position_col, char_col=control_col)
        properties = np.zeros((len(property_sets),len(self.WT_seq)))
        properties[:,:]=np.nan
        p_tick_pos = []
        p_tick_lab = []
        for i, key in enumerate(property_sets.keys()):
            p_tick_pos.append(i)
            p_tick_lab.append(key[0])
            to_set = to_collapse.loc[key].reset_index()
            properties[np.ones(to_set.shape[0]).astype(int)*i,to_set[position_col].values.astype(int)-1] = to_set['Statistic'].values
        
        if norm is None:
            norm=CenteredNorm(vcenter=0)

        if grid:
            plt.pcolormesh(properties,shading='flat', edgecolors='k',
                    cmap=cmap, linewidth=grid, norm=norm)
            
            ax.invert_yaxis()
            offset=0.5
            # axes[0,0].set_aspect('equal')
        else:
            plt.imshow(properties,interpolation='none', 
                    cmap=cmap, norm=norm)
            offset=0
        ax.xaxis.set_visible(False)
        cbar = plt.colorbar(orientation='vertical',cax=cax)
        cax.set_yscale('linear')
        ax.set_yticks(np.array(p_tick_pos)+offset,p_tick_lab)
        ax.set_xlim(plot_pos)

    def draw_feature_sets(self, ax,feature_sets, grid, plot_pos, seq_start, allow_multi_rows=False):
        """Draw heatmap of feature sets.
        
        Args:
            ax (matplotlib.axes.Axes): Axes to draw heatmap on
            feature_sets (list): List of dictionaries containing feature information
            grid (bool): Whether to draw grid lines
            plot_pos (tuple): (start, end) positions to plot
            seq_start (int): Starting position in sequence
            allow_multi_rows (bool, optional): Whether to allow multiple rows per feature. Defaults to False.
        """
        plt.sca(ax)
        if allow_multi_rows:
            num_features = len(np.unique([feature['name'] for feature in feature_sets]))
        else:
            num_features = len(feature_sets)
        feature_map = np.zeros((num_features,len(self.WT_seq)))

        cmap = ['xkcd:white']
        names = []
        i=0
        for j, feature in enumerate(feature_sets):
            feature_pos = np.array(feature['positions'])-seq_start
            assert (feature_pos>=0).all(), 'Features included with positions before the start of the sequence'
            if allow_multi_rows and (feature['name'] in names):
                feature_map[names.index(feature['name']), feature_pos]=j+1
            else:
                feature_map[i, feature_pos]=j+1
                i+=1
                names.append(feature['name'])
            cmap.append(feature['color'])
            

        cmap = sns.blend_palette(cmap,as_cmap=True)
        if grid:
            plt.pcolormesh(feature_map,shading='flat', edgecolors='k',
                    cmap=cmap,linewidth=grid)
            
            ax.invert_yaxis()
            offset=0.5
            # axes[1,0].grid(True, c='k')
        else:
            plt.imshow(feature_map,interpolation='none', aspect='equal',
                    cmap=cmap,)
            offset=0

        ax.xaxis.set_visible(False)
        ax.set_yticks(np.arange(len(names))+offset,names)
        ax.yaxis.tick_right()
        ax.set_xlim(plot_pos)

    
    def make_singles_heatmap(self, mut_num_col='Mut_num', 
                             position_col='Position', main_cmap='RdBu', mutability_cmap=sns.blend_palette(['#CD3333','xkcd:white'],as_cmap=True),
                               bad_color='xkcd:grey', sns_context='notebook', rc_params={}, grid=False,invert_mutability=False, norm = None, seq_start=0,
                               tick_interval=5, plot_pos = None, plot_mutability=False, plot_properties=False, plot_features=False, mutability_norm=None,
                               property_norm=None,
                               property_sets = {'Hydrophobic':['A','V','I','L','M','F','Y','W'],'Polar':['S','T','N','Q'],'Charged':['R','H','K','D', 'E']},
                               property_cmap='RdBu',
                               feature_sets=[],
                               control_col='AA', WT_vals=['WT'], stop_vals=['*'], thresh=0,allow_multi_rows=False
                               ):
        """Create a comprehensive heatmap visualization of single mutant effects.
        
        Args:
            mut_num_col (str, optional): Column containing number of mutations. Defaults to 'Mut_num'.
            position_col (str, optional): Column containing positions. Defaults to 'Position'.
            main_cmap (str or matplotlib.colors.Colormap, optional): Colormap for main heatmap. Defaults to 'RdBu'.
            mutability_cmap (matplotlib.colors.Colormap, optional): Colormap for mutability scores.
            bad_color (str, optional): Color for missing values. Defaults to 'xkcd:grey'.
            sns_context (str, optional): Seaborn context. Defaults to 'notebook'.
            rc_params (dict, optional): RC parameters. Defaults to {}.
            grid (bool, optional): Whether to show grid lines. Defaults to False.
            invert_mutability (bool, optional): Whether to invert mutability scores. Defaults to False.
            norm (matplotlib.colors.Normalize, optional): Normalization for main heatmap. Defaults to None.
            seq_start (int, optional): Starting position in sequence. Defaults to 0.
            tick_interval (int, optional): Interval between ticks. Defaults to 5.
            plot_pos (tuple, optional): (start, end) positions to plot. Defaults to None.
            plot_mutability (bool, optional): Whether to plot mutability scores. Defaults to False.
            plot_properties (bool, optional): Whether to plot property enrichments. Defaults to False.
            plot_features (bool, optional): Whether to plot feature sets. Defaults to False.
            mutability_norm (matplotlib.colors.Normalize, optional): Normalization for mutability scores. Defaults to None.
            property_norm (matplotlib.colors.Normalize, optional): Normalization for property scores. Defaults to None.
            property_sets (dict, optional): Dictionary of property sets. Defaults to standard amino acid properties.
            property_cmap (str or matplotlib.colors.Colormap, optional): Colormap for properties. Defaults to 'RdBu'.
            feature_sets (list, optional): List of feature dictionaries. Defaults to [].
            control_col (str, optional): Column containing controls. Defaults to 'AA'.
            WT_vals (list, optional): Values indicating wild type. Defaults to ['WT'].
            stop_vals (list, optional): Values indicating stop codons. Defaults to ['*'].
            thresh (float, optional): Threshold for calculations. Defaults to 0.
            allow_multi_rows (bool, optional): Allow multiple rows per feature. Defaults to False.
            
        Returns:
            tuple: (matplotlib.figure.Figure, dict of matplotlib.axes.Axes)
        """
        if plot_pos is None:
            plot_pos=(0, len(self.WT_seq))
        else:
            plot_pos = (plot_pos[0]-seq_start, plot_pos[1]-seq_start)
            
        with rc_context(rc_params), sns.plotting_context(sns_context):
            

            #Setup figure
            heights = []
            panel_spec = []

            if plot_mutability:
                heights.append(1)
                panel_spec.append(['mutability_hmap','mutability_cmap'])

            if plot_properties:
                heights.append(len(property_sets))
                panel_spec.append(['property_hmap','property_cmap'])

            if plot_features:
                heights.append(len(feature_sets))
                panel_spec.append(['feature_hmap','blank1'])

            heights.append(len(self.alphabet)+1)
            panel_spec.append(['mutant_hmap','mutant_cmap'])
            heights.append(0.5)
            panel_spec.append(['buffer','buffer'])
            heights.append(1)
            panel_spec.append(['mutant_second_axis','blank'])

            width_ratios=[plot_pos[-1],plot_pos[-1]//50]
            gs_kw = dict(width_ratios=width_ratios, height_ratios=heights)
            fig, axd = plt.subplot_mosaic(panel_spec,
                                        gridspec_kw=gs_kw)



            #Plot Panels
            axd['blank'].patch.set_visible(False)
            axd['blank'].yaxis.set_visible(False)
            axd['blank'].xaxis.set_visible(False)
            axd['blank'].spines[['right', 'top', 'left', 'bottom']].set_visible(False)
            axd['buffer'].patch.set_visible(False)
            axd['buffer'].yaxis.set_visible(False)
            axd['buffer'].xaxis.set_visible(False)
            axd['buffer'].spines[['right', 'top', 'left', 'bottom']].set_visible(False)

            if isinstance(main_cmap, str):
                from matplotlib import colormaps
                main_cmap = colormaps.get_cmap(main_cmap)
            else:
                main_cmap = main_cmap
            main_cmap.set_bad(color=bad_color)

            self.draw_singles_heatmap(axd['mutant_hmap'], axd['mutant_cmap'], axd['mutant_second_axis'],
                                       mut_num_col, control_col, WT_vals, stop_vals, 
                                      position_col, grid, main_cmap, norm,
                                      seq_start, tick_interval, plot_pos)
            
            main_pos = axd['mutant_hmap'].get_position().bounds

            if plot_mutability:
                if isinstance(mutability_cmap, str):
                    from matplotlib import colormaps
                    mutability_cmap = colormaps.get_cmap(mutability_cmap)
                else:
                    mutability_cmap = mutability_cmap
                mutability_cmap.set_bad(color=bad_color)
                self.draw_mutability_score(axd['mutability_hmap'], axd['mutability_cmap'], 
                                           control_col, WT_vals, stop_vals, position_col, 
                                           thresh,invert_mutability, grid, mutability_cmap, mutability_norm,
                                           plot_pos)
                new_pos = axd['mutability_hmap'].get_position().bounds
                axd['mutability_hmap'].set_position((main_pos[0],new_pos[1],main_pos[2],new_pos[3]))

            if plot_properties:
                if isinstance(property_cmap, str):
                    from matplotlib import colormaps
                    property_cmap = colormaps.get_cmap(property_cmap)
                else:
                    property_cmap = property_cmap
                property_cmap.set_bad(color=bad_color)

                self.draw_property_scores(axd['property_hmap'], axd['property_cmap'],
                                          property_sets, thresh, position_col, control_col, 
                                          property_cmap, grid, property_norm,
                                          plot_pos)
                new_pos = axd['property_hmap'].get_position().bounds
                axd['property_hmap'].set_position((main_pos[0],new_pos[1],main_pos[2],new_pos[3]))

            if plot_features:
                self.draw_feature_sets(axd['feature_hmap'],feature_sets, grid, plot_pos,seq_start,allow_multi_rows=allow_multi_rows)
                axd['blank1'].patch.set_visible(False)
                axd['blank1'].yaxis.set_visible(False)
                axd['blank1'].xaxis.set_visible(False)
                axd['blank1'].spines[['right', 'top', 'left', 'bottom']].set_visible(False)
                new_pos = axd['feature_hmap'].get_position().bounds
                axd['feature_hmap'].set_position((main_pos[0],new_pos[1],main_pos[2],new_pos[3]))
            
            
            return(fig,axd)

