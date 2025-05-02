from matplotlib import rc_context
import seaborn as sns
from typing import Union
rc_defaults = {
   "font.family" : 'sans-serif',
   "font.sans-serif" : 'Helvetica',
   'figure.autolayout' :  True,
   "figure.figsize" : (1.4,1.4),
   "figure.dpi": 400
}

sns_defaults = {
    'axes.linewidth' : 0.25,
    'grid.linewidth' : 0.25,
    'lines.linewidth' : 0.5,
    'lines.markersize' : 1,
    'patch.linewidth' : 0.25,
    'xtick.major.width' : 0.25,
    'ytick.major.width' : 0.25,
    'xtick.minor.width' : 0.25,
    'ytick.minor.width' : 0.25,
    'xtick.major.size' : 1,
    'ytick.major.size' : 1,
    'xtick.minor.size' : 0.75,
    'xtick.major.width' : 0.25,
    'ytick.major.width' : 0.25,
    'xtick.minor.width' : 0.25,
    'ytick.minor.width' : 0.25,
    'xtick.major.size' : 1,
    'ytick.major.size' : 1,
    'xtick.minor.size' : 0.75,
    'ytick.minor.size' : 0.75,
    'font.size' : 6,
    'axes.labelsize' : 6,
    'axes.titlesize' : 6,
    'xtick.labelsize' : 6,
    'ytick.labelsize' : 6,
    'legend.fontsize' : 6,
    'legend.title_fontsize' : 6
}



class PlottingContext:
    """Context manager for setting matplotlib and seaborn plotting parameters.
    
    Args:
        context (str): The plotting context - one of 'paper', 'notebook', 'talk', 'poster'
        font_scale (float): Scale factor for font sizes
        rc (dict): Dictionary of rc parameters to override
    """
    def __init__(self, sns_settings: Union[dict, str] = None, rc_settings: dict = None, palette: Union[str, list, sns.palettes._ColorPalette] = None):
        # Start with defaults and update with any passed values
        self.sns_settings = sns_defaults.copy()
        if sns_settings:
            if isinstance(sns_settings, str):
                self.sns_settings = sns_settings
            else:
                self.sns_settings.update(sns_settings)

        self.sns_settings = sns.plotting_context(self.sns_settings)
            
        self.rc_settings = rc_defaults.copy() 
        if rc_settings:
            self.rc_settings.update(rc_settings)

        self.rc_settings = rc_context(self.rc_settings)
        
        # Store the original palette
        if not isinstance(palette, sns.palettes._ColorPalette):
            self.palette = sns.color_palette(palette)
        else:
            self.palette = palette
        self.original_palette = None
    def __enter__(self):
        self.rc_context = rc_context(self.rc_settings)
        self.sns_settings.__enter__()
        self.rc_settings.__enter__()
        
        if self.palette is not None:
            # Store the original palette
            self.original_palette = sns.color_palette()
            sns.set_palette(self.palette)
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        # Reset to original palette
        if self.original_palette is not None:
            sns.set_palette(self.original_palette)
            self.original_palette = None
        self.rc_settings.__exit__(exc_type, exc_val, exc_tb)
        self.sns_settings.__exit__(exc_type, exc_val, exc_tb)
