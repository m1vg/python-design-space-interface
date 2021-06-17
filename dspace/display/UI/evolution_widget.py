import dspace
import pandas as pd
from decimal import Decimal
from math import log10, floor
from dspace.SWIG.dspace_interface import *
import numpy as np
from distutils.version import LooseVersion, StrictVersion
import IPython
np.seterr(all='raise')
from functools import partial
import matplotlib.pyplot as plt
import cStringIO
from matplotlib.backends.backend_agg import FigureCanvasAgg
import itertools


if StrictVersion(IPython.__version__) < StrictVersion('4.0.0'):
    from IPython.html.widgets import interact, interactive, fixed
    from IPython.html.widgets import HTMLWidget as HTML
    from IPython.html.widgets import TabWidget as Tab
    from IPython.html.widgets import CheckboxWidget as Checkbox
    from IPython.html.widgets import ButtonWidget as Button
    from IPython.html.widgets import ContainerWidget as Box
    from IPython.html.widgets import TextWidget as Text
    from IPython.html.widgets import TextareaWidget as Textarea
    from IPython.html.widgets import DropdownWidget as Dropdown
    from IPython.html.widgets import RadioButtonsWidget as RadioButtons
    from IPython.html.widgets import PopupWidget as Popup
    from IPython.html.widgets import LatexWidget as Latex
    from IPython.html.widgets import FloatTextWidget as FloatText
    from IPython.html.widgets import ImageWidget as Image
    from IPython.html.widgets import IntSliderWidget as Slider
    from IPython.html.widgets import FloatSliderWidget as FloatSlider
    VBox = Box
    HBox = Box
    old_ipython = True
else:
    from ipywidgets import *
    from popup import Popup
    old_ipython = False


class Evolution(object):

    def __init__(self, controller):
        setattr(self, 'controller', controller)

    def update_task_fields(self, name, value):

        if value == 'Plots':
            self.phenotypic_plots.visible = True
            self.phenotypic_neighbors.visible = False
        else:
            self.phenotypic_plots.visible = False
            self.phenotypic_neighbors.visible = True

    def create_evolution_widget(self):

        controller = self.controller
        if controller.ds is None:
            return

        self.phenotypic_neighbors = Box()
        self.phenotypic_plots = Box()
        self.phenotypic_plots.visible = False

        task = Dropdown(description='Task',
                        values=['Mutation Rates/Neighbors', 'Plots'],
                        options=['Mutation Rates/Neighbors', 'Plots'],
                        value='Mutation Rates/Neighbors')

        task.on_trait_change(self.update_task_fields, 'value')

        wi = Box(
                children=[task,
                          self.phenotypic_neighbors,
                          self.phenotypic_plots
                          ]
                 )

        self.update_phenotypic_neighbors()
        self.update_phenotypic_plots()

        return ('Evolutionary Studies', wi)

    def update_phenotypic_plots(self):

        title = HTML(value='<b>Plots  </b>')
        plot_type = Dropdown(description='Plot Type',
                             values=['Phenotype Frequencies'],
                             options=['Phenotype Frequencies'],
                             value='Phenotype Frequencies')

        phenotypes = Text(description='Phenotype List:',
                          value='1, 3, 5, 6, 7, 8, 11, 15, 16')

        valid_phenotypes = Checkbox(description='Valid Phenotypes', value=False)
        fun = lambda name, value: self.update_phenotype_field(name, value, phenotypes=phenotypes)
        valid_phenotypes.on_trait_change(fun, 'value')
        valid_phenotypes.value = True

        method = Dropdown(description='Method',
                          values=['Tolerances',
                                  'Bounding Box',
                                  'Geometric Mean T. & BB.',
                                  'Vertex Enumeration'],
                          options=['Tolerances',
                                  'Bounding Box',
                                  'Geometric Mean T. & BB.',
                                  'Vertex Enumeration'],
                          value='Tolerances',
                          )

        color = Dropdown(description='Color',
                         values=['black', 'red', 'darkgreen', 'darkcyan', 'gold', 'grey'],
                         options=['black', 'red', 'darkgreen', 'darkcyan', 'gold', 'grey'],
                         value='black')

        methods_box = Box(children=[method,
                                    color])

        self.methods_list = [{'method': method, 'color': color}]
        add_button = Button(description='Add Method')
        delete_button = Button(description='Delete Method', visible=False)
        delete_button.methods_box = methods_box

        lower_bound = Slider(description='Lower Bound Par. Space (Log)',
                               value=-3,
                               min=-9,
                               max=-1,
                               step=1,
                               visible=True)

        upper_bound = Slider(description='Upper Bound Par. Space (Log)',
                               value=3,
                               min=1,
                               max=9,
                               step=1,
                               visible=True)

        bounds_box = Box(children=[lower_bound,
                                   upper_bound
                                  ])

        add_button.methods_box = methods_box
        add_button.delete_button = delete_button
        add_button.on_click(self.update_methods_list)
        delete_button.on_click(self.update_methods_list)

        add_plot = Button(description='Add Plot')
        add_plot.lower_bound = lower_bound
        add_plot.upper_bound = upper_bound
        add_plot.phenotypes = phenotypes
        add_plot.valid_phenotypes = valid_phenotypes

        add_plot.on_click(self.make_phenotype_frequencies_plot)

        controls_box = Box(children=[delete_button,
                                     add_button,
                                     ])

        steady_state_distribution = Box()
        header_plots = Box(children=[title,
                                     plot_type])

        steady_state_distribution.children = [phenotypes,
                                              valid_phenotypes,
                                              methods_box,
                                              controls_box,
                                              bounds_box,
                                              add_plot]

        self.phenotypic_plots.steady_state_distribution = steady_state_distribution
        self.phenotypic_plots.header_plots = header_plots
        self.update_response_times_plot()
        plot_type.on_trait_change(self.update_plot_menu, 'value')

    def update_plot_menu(self, name, value):
        if value == 'Phenotype Frequencies':
            self.phenotypic_plots.steady_state_distribution.visible = True
            self.phenotypic_plots.response_time.visible = False
        else:
            self.phenotypic_plots.steady_state_distribution.visible = False
            self.phenotypic_plots.response_time.visible = True

    def update_response_times_plot(self):
        controller = self.controller

        # Get the auxiliary, independent and dependent variables
        auxiliary = controller.ds.auxiliary_variables
        dependent = controller.ds.dependent_variables
        dependent_no_artificial = controller.ds.dependent_variables

        # 2. Remove auxiliary variables from dependent variables
        for aux in auxiliary:
            if aux in dependent:
                dependent.remove(aux)

        # Generate conserved, artificial variables Xc1, ... Xcn
        if controller.ds.number_conservations != 0:
            artificial_var = ['Xc' + str(n+1) for n in range(int(controller.ds.number_conservations))]

            # Remove artificial variables from dependent variables
            for aux in artificial_var:
                if aux in dependent_no_artificial:
                    dependent_no_artificial.remove(aux)

        response_time_box = Box()

        x_axis = Dropdown(description='X-Axis')
        x_axis.values = dict((k, k) for k in controller.ds.independent_variables)

        x_min = FloatText(description='X-Min', value=0.01, visible=False)
        x_max = FloatText(description='X-Max', value=100, visible=False)

        custom_points_x_axis = Checkbox(value=False, description='Custom Values for X-axis')
        custom_points_titration = Textarea(description='Values',
                                           value='1, 2, 3',
                                           visible=True)

        y_axis = Dropdown(description="Y-Axis")
        y_axis.values = dict((k, k) for k in dependent_no_artificial)

        resp_time_method = Dropdown(description='Response Time Method')
        resp_time_method.values = dict((k, k) for k in ['y_final - y_initial', 'max(y) - min(y)', 'custom'])

        reference_state = Dropdown(description='Reference State')
        reference_state.values = dict((k, k) for k in ['First Point', 'Previous Point'])

        initial_string = HTML(value='<b>Initial Concentrations</b>')
        self.initial_concentration = dict((i, FloatText(description=i, value=0.001)) for i in dependent)
        add_plot = Button(description='Add Plot')

        response_time_box.children = [x_axis,
                                      x_min,
                                      x_max,
                                      custom_points_x_axis,
                                      custom_points_titration,
                                      y_axis,
                                      resp_time_method,

                                      reference_state,
                                      initial_string] + \
                                      [value for key, value in self.initial_concentration.iteritems()] + \
                                      [add_plot
                                      ]

        self.phenotypic_plots.response_time = response_time_box

        self.phenotypic_plots.children = [  self.phenotypic_plots.header_plots,
                                            self.phenotypic_plots.steady_state_distribution,
                                            self.phenotypic_plots.response_time
                                          ]

        response_time_box.visible = False

    def make_phenotype_frequencies_plot(self, b):

        # 1. Determine the list of phenotypes for which the plot will be generated.

        if b.valid_phenotypes.value is True:
            phenotypes = self.controller.ds.valid_cases()
        else:
            phenotypes = b.phenotypes.split(',')

        # Initialice plots
        fig = plt.figure(figsize=[5.5698, 4], dpi=600, facecolor='w')
        ax = fig.add_axes([0.1714, 0.2, 0.6, 0.7])
        ax.set_title('Phenotype Frequencies Plot')
        ax.set_xlabel('Phenotype Number, k')

        # for each method, generate plot
        for method in self.methods_list:
            self.make_plot_for_method(ax, b, method=method, phenotypes=phenotypes)

        # 3. Make the plots

        buf = cStringIO.StringIO()
        canvas = FigureCanvasAgg(fig)
        canvas.print_png(buf)
        data = buf.getvalue()
        plt.close()
        self.controller.figures.add_figure(data,
                                      title='Analysis of the ' + self.controller.name + ' by system design space showing a phenotype frequencies plot',
                                      caption='Nominal distribution of phenotype frequencies based only on analytically determined volumes.',
                                      pvals=None,
                                      colors=None)

    def make_plot_for_method(self, ax, b, method=None, phenotypes=None):

        lb = 10 ** b.lower_bound.value
        ub = 10 ** b.upper_bound.value
        volume_method = method['method'].value
        color = method['color'].value

        # Determine normalized volume for each phenotype
        total_vol = 0
        vol_list = []
        for c in phenotypes:
            phenotype = self.controller.ds(str(c))
            volume = phenotype.volume(
                                            ignore_unbounded=False,
                                            log_coordinate=True,
                                            shared_boundaries=False,
                                            method=volume_method,
                                            lowerBounds=lb,
                                            upperBounds=ub,
                                            maxVertices=0,
                                            limitVertices=False
                                            )
            total_vol += volume
            vol_list.append(volume)

        norm_vol = []
        for vol in vol_list:
            norm_vol.append(vol/total_vol)

        time = np.linspace(1, len(phenotypes), len(phenotypes))
        ax.plot(time, norm_vol, marker='o', color=color)
        ax.set_ylabel('Frequency, Rk')
        ax.set_xticks(time)
        ax.set_xticklabels(phenotypes)
        ax.set_yscale("log")

    def update_phenotype_field(self, name, value, phenotypes=None):
        if value is True:
            if phenotypes is not None:
                phenotypes.visible = False
        else:
            if phenotypes is not None:
                phenotypes.visible = True

    def update_methods_list(self, b):

        if b.description == 'Add Method':

            method = Dropdown(description='Method',
                              values=['Tolerances',
                                      'Bounding Box',
                                      'Geometric Mean T. & BB.',
                                      'Vertex Enumeration'],
                              options=['Tolerances',
                                       'Bounding Box',
                                       'Geometric Mean T. & BB.',
                                       'Vertex Enumeration'],
                              value='Tolerances',
                              )

            color = Dropdown(description='Color',
                             values=['black', 'red', 'darkgreen', 'darkcyan', 'gold', 'grey'],
                             options=['black', 'red', 'darkgreen', 'darkcyan', 'gold', 'grey'],
                             value='black')

            self.methods_list.append({'method': method, 'color': color})

            children = []
            for element in self.methods_list:
                children.append(element['method'])
                children.append(element['color'])
            b.methods_box.children = children
            b.delete_button.visible = True

        elif b.description == 'Delete Method':

            self.methods_list.pop()
            children = []
            for element in self.methods_list:
                children.append(element['method'])
                children.append(element['color'])
                children.append(HTML(value='</br>'))
            b.methods_box.children = children
            if len(self.methods_list) == 1:
                b.visible = False

    def update_phenotypic_neighbors(self, bol_intersecting=False, bol_share=True, phen_value='', method='Signature'):

        title = HTML(value='<b> Phenotypic Neighbors & Phenotype-specific Mutation Rates  </b>')
        phenotype = Text(description='* Calculate Neighbors / Mutation Rates for Phenotype Number: ', value=phen_value)
        intersecting = Checkbox(description="Intersects with Phenotype",
                                value=bol_intersecting,
                                visible=False)
        explicitly_share = Checkbox(description="Has a Common Boundary with Phenotype",
                                    value=bol_share,
                                    visible=False)
        export_neighbors = Checkbox(description='Export Neighbors (.xlsx)', value=False)
        neighbors_method = Dropdown(description='Method',
                                    values=['Vertex Enumeration',
                                            'Signature',
                                            'User Defined',
                                            'All Mutation Rates'],
                                    options=['Vertex Enumeration',
                                             'Signature',
                                             'User Defined',
                                             'All Mutation Rates'],
                                    value=method
                                    )

        user_cases = Textarea(description="User Specified:",
                                  visible=False)

        neighbors_method.on_trait_change(self.update_column_block_neighbors, 'value')

        nr_mutations = Dropdown(description='# Mutations Allowed',
                                 values=['Multiple',
                                         'Single'],
                                 options=['Multiple',
                                          'Single'],
                                 value='Multiple',
                                 visible=True
                                 )

        extra_columns = VBox()
        self.extra_columns = extra_columns
        filters = VBox(children=[])

        all_mutation_rates_block_dic = {'nr_mutations': nr_mutations,                                                                                   # 0

                                        'operating_point_calculation_method': Dropdown(description='Operating Point Calculation Method',                # 1
                                                                                         values=[
                                                                                                 'Tolerances',
                                                                                                 'Bounding Box',
                                                                                                 'Geometric Mean T. & BB.',
                                                                                                 'Vertex Enumeration',
                                                                                                ],
                                                                                         options=[
                                                                                                  'Tolerances',
                                                                                                  'Bounding Box',
                                                                                                  'Geometric Mean T. & BB.',
                                                                                                  'Vertex Enumeration',
                                                                                                    ],
                                                                                         value='Tolerances',
                                                                                         visible=True,
                                                                                         ),

                                        'volume_calculation_method': Dropdown(description='Volume Calculation Method',                                  # 2
                                                                                 values=[
                                                                                         'Tolerances',
                                                                                         'Bounding Box',
                                                                                         'Geometric Mean T. & BB.',
                                                                                         'Vertex Enumeration',
                                                                                         ],
                                                                                 options=[
                                                                                         'Tolerances',
                                                                                         'Bounding Box',
                                                                                         'Geometric Mean T. & BB.',
                                                                                         'Vertex Enumeration',
                                                                                          ],
                                                                                 value='Tolerances',
                                                                                 visible=True
                                                                                 ),

                                        'tracks_calculation_method': Dropdown(description='Tracks Calculation Method',                                  # 3
                                                                             values=['Vertex Enumeration',
                                                                                     'Bounding Box'],
                                                                             options=['Vertex Enumeration',
                                                                                      'Bounding Box'],
                                                                             value='Bounding Box',
                                                                             visible=False
                                                                             ),

                                        'method_for_width_determination': Dropdown(description='Method for Width Determination',                        # 4
                                                                                     values=['Width at Orthogonal Centroid',
                                                                                             'Average Width',
                                                                                             'Average of Methods',
                                                                                             'Average Width over Facets at Centroid',
                                                                                             'Custom Formulas'],

                                                                                     options=['Width at Orthogonal Centroid',
                                                                                              'Average Width',
                                                                                              'Average of Methods',
                                                                                              'Average Width over Facets at Centroid',
                                                                                              'Custom Formulas'],
                                                                                     value='Average Width over Facets at Centroid',
                                                                                     visible=False),

                                        'enforce_bounds_operating_point': Checkbox(description='Enforce Bounds Operating Point',
                                                                                     value=True,
                                                                                     visible=True),                                                    # 5

                                        'enforce_bounds_volume': Checkbox(description='Enforce Bounds Volume',
                                                                             value=True,
                                                                             visible=True),                                                            # 6

                                        'lower_bound_par_space': Slider(description='Lower Bound Par. Space (Log)',                                    # 7
                                                                       value=-3,
                                                                       min=-9,
                                                                       max=-1,
                                                                       step=1,
                                                                       visible=True),

                                        'upper_bound_par_space': Slider(description='Upper Bound Par. Space (Log)',                                    # 8
                                                                           value=3,
                                                                           min=1,
                                                                           max=9,
                                                                           step=1,
                                                                           visible=True),

                                        'delta_value': FloatText(description='Delta Value',                                                             # 9
                                                                  value=2**0.5,
                                                                  visible=True),

                                        'lambda_value': FloatText(description='Lambda Value',                                                           # 10
                                                                  value=2.0,
                                                                  visible=True),

                                        'include_fractional_volume_data': Checkbox(description='Include Fractional Volume Data',                        # 11
                                                                                     value=False,
                                                                                     visible=False),

                                        'include_real_volumes': Checkbox(description='Include Real Volumes (lrs)',                                      # 12
                                                                             value=False,
                                                                             visible=False),

                                        'file_name': Text(description='File Name ',                                                                     # 13
                                                          value=self.controller.name + '_Parameters.xlsx'),

                                        'asymmetrical_bounds_checkbox': Checkbox(description='Asymmetrical Bounds', value=False),
                                        'asymmetrical_bounds_box': VBox(visible=False),
                                        'asymmetrical_bounds_list': [],
                                        }

        all_mutation_rates_block = VBox(children=[all_mutation_rates_block_dic['nr_mutations'],                                             # 0
                                                  all_mutation_rates_block_dic['operating_point_calculation_method'],                       # 1
                                                  all_mutation_rates_block_dic['volume_calculation_method'],                                # 2
                                                  all_mutation_rates_block_dic['tracks_calculation_method'],                                # 3
                                                  all_mutation_rates_block_dic['method_for_width_determination'],                           # 4
                                                  all_mutation_rates_block_dic['enforce_bounds_operating_point'],                           # 5
                                                  all_mutation_rates_block_dic['enforce_bounds_volume'],                                    # 6


                                                  all_mutation_rates_block_dic['asymmetrical_bounds_checkbox'],
                                                  all_mutation_rates_block_dic['asymmetrical_bounds_box'],


                                                  all_mutation_rates_block_dic['lower_bound_par_space'],                                    # 7
                                                  all_mutation_rates_block_dic['upper_bound_par_space'],                                    # 8
                                                  all_mutation_rates_block_dic['delta_value'],                                              # 9
                                                  all_mutation_rates_block_dic['lambda_value'],                                             # 10
                                                  all_mutation_rates_block_dic['include_fractional_volume_data'],                           # 11
                                                  all_mutation_rates_block_dic['include_real_volumes'],                                     # 12
                                                  all_mutation_rates_block_dic['file_name']                                                 # 13

                                                 ],

                                        visible=False,
                                        )

        nr_mutations.on_trait_change(self.update_all_mutations_block, 'value')

        el = ['volume_calculation_method', 'include_fractional_volume_data', 'asymmetrical_bounds_checkbox']
        for i in el:
            all_mutation_rates_block_dic[i].on_trait_change(self.update_all_mutations_block, 'value')

        b = Button(description='Calculate')
        b.title = title
        b.phenotype = phenotype
        b.explicitly_share = explicitly_share
        b.intersecting = intersecting
        b.on_click(self.get_neighbors)
        b.extra_columns = extra_columns
        b.filters = filters
        b.export_table = export_neighbors
        b.neighbors_method = neighbors_method
        b.user_cases = user_cases
        b.all_mutation_rates = all_mutation_rates_block

        remove_column = Button(description='Remove column')
        remove_column.visible = False
        remove_column.on_click(self.remove_neighbors_column)
        remove_column.column_block = extra_columns

        add_column = Button(description='Add Column')
        add_column.on_click(self.add_neighbors_column)
        add_column.column_block = extra_columns
        add_column.remove_column = remove_column

        add_filter = Button(description='Add Column Filter')
        add_filter.on_click(self.add_neighbors_filter)
        add_filter.column_block = extra_columns
        add_filter.filters = filters

        phenotypic_neighbors_block = VBox()
        phenotypic_neighbors_block.children = [title,
                                               phenotype,
                                               neighbors_method,
                                               intersecting,
                                               explicitly_share,
                                               user_cases,
                                               all_mutation_rates_block,

                                               add_column,
                                               extra_columns,
                                               remove_column,

                                               filters,
                                               add_filter,
                                               export_neighbors,

                                               b]

        self.phenotypic_neighbors.children = [phenotypic_neighbors_block
                                              ]

        self.phenotypic_neighbors.title = title
        self.phenotypic_neighbors.phenotype = phenotype
        self.phenotypic_neighbors.explicitly_share = explicitly_share
        self.phenotypic_neighbors.intersecting = intersecting
        self.phenotypic_neighbors.export_table = export_neighbors
        self.phenotypic_neighbors.neighbors_method = neighbors_method
        self.phenotypic_neighbors.user_cases = user_cases
        self.phenotypic_neighbors.all_mutation_rates = all_mutation_rates_block_dic
        self.phenotypic_neighbors.all_mutation_rates_vbox = all_mutation_rates_block
        self.populate_asymmetrical_box()
        self.phenotypic_neighbors.add_column = add_column
        self.phenotypic_neighbors.add_filter = add_filter
        self.phenotypic_neighbors.action_button = b

    def populate_asymmetrical_box(self):
        block = self.phenotypic_neighbors.all_mutation_rates

        ind_variables = self.controller.ds.independent_variables
        asymm_list = []
        children = []

        for variable in ind_variables:
            asymm_dic = {}
            asymm_dic['name'] = variable
            asymm_dic['low_bound'] = FloatSlider(description='Lower Bound (Log)',
                                                 value=-3,
                                                 min=-9,
                                                 max=0,
                                                 step=0.25)

            asymm_dic['upper_bound'] = FloatSlider(description='Upper Bound (Log)',
                                                   value=3,
                                                   min=0,
                                                   max=9,
                                                   step=0.25,
                                                   )
            asymm_list.append(asymm_dic)
            children.append(HTML(value='<b>Bounds for ' + asymm_dic['name'] + '</b>'))
            children.append(asymm_dic['low_bound'])
            children.append(asymm_dic['upper_bound'])

        block['asymmetrical_bounds_list'] = asymm_list
        block['asymmetrical_bounds_box'].children = children


    def remove_neighbors_column(self, b):
        controller = self.controller
        if len(b.column_block.children) == 0:
            return
        if len(b.column_block.children) == 1:
            b.visible = False
        children = [i for i in b.column_block.children[:-1]]
        b.column_block.children = children

    def add_neighbors_column(self, b):
        controller = self.controller
        children = [i for i in b.column_block.children]

        nr_ineq = DSDesignSpaceNumberOfBoundaries(controller.ds._swigwrapper)*2 + len(controller.constraints) + \
                  len(controller.pvals)*2
        max_vertices = nr_ineq ** (len(controller.pvals) / 2)

        new_column = [Dropdown(description='Column ' + str(len(children) + 3),                          # 0
                               values=[
                                        'Distance',
                                        'Number of Vertices',
                                        'Dimension of Shared Space',
                                        'Number of Shared Boundaries',
                                        'Volume of Shared Space',
                                        'Intersect',
                                        'Mutation Rate'],

                               options=[
                                        'Distance',
                                        'Number of Vertices',
                                        'Dimension of Shared Space',
                                        'Number of Shared Boundaries',
                                        'Volume of Shared Space',
                                        'Intersect',
                                        'Mutation Rate',
                                        ],

                               value='Distance'),

                      Dropdown(description='Volume Calculation Method',                                 # 1
                               values=[
                                        'Tolerances',
                                        'Bounding Box',
                                        'Geometric Mean T. & BB.',
                                        'Vertex Enumeration'
                                ],

                               options=[
                                        'Tolerances',
                                        'Bounding Box',
                                        'Geometric Mean T. & BB.',
                                        'Vertex Enumeration',
                               ],
                               value='Tolerances',
                               visible=False
                               ),

                      Checkbox(description='Limit Number of Vertices',                                  # 2
                               value=False,
                               visible=False),

                      Slider(description='Max. Number of Vertices',                                     # 3
                             value=max_vertices,
                             min=1,
                             max=round_sig(max_vertices, 1)*2,
                             step=round_sig(max_vertices, 1)/100,
                             visible=False
                             ),

                      Dropdown(description='Operating Point Calculation Method',                        # 4
                               values=[
                                       'Tolerances',
                                       'Bounding Box',
                                       'Geometric Mean T. & BB.',
                                       'Vertex Enumeration',
                                       'Grid',
                               ],
                               options=[
                                       'Tolerances',
                                       'Bounding Box',
                                       'Geometric Mean T. & BB.',
                                       'Vertex Enumeration',
                                       'Grid',
                               ],
                               value='Tolerances',
                               visible=False,
                               ),

                      Checkbox(description='Enforce Bounds', value=True, visible=False),               # 5

                      Checkbox(description='Asymmetrical Bounds', value=False, visible=False),          # 6

                      Box(children=[                                                                    # 7
                                      Slider(description='Lower Bound Par. Space (Log)',                # 7.0
                                             value=-3,
                                             min=-9,
                                             max=-1,
                                             step=1,
                                             visible=False),

                                      Slider(description='Upper Bound Par. Space (Log)',               # 7.1
                                             value=3,
                                             min=1,
                                             max=9,
                                             step=1,
                                             visible=False),
                      ]),

                      Dropdown(description='Ignore Unbounded Parameters',                               # 8
                               values=[True, False],
                               options=[True, False],
                               value=False,
                               visible=False),

                      Dropdown(description='Logarithmic Coordinates',                                   # 9
                               values=[True, False],
                               options=[True, False],
                               value=True,
                               visible=False),

                      FloatText(description='Delta Value',                                               # 10
                                value=2**0.5,
                                visible=False),

                      FloatText(description='Lambda Value',                                             # 11
                                value=2.0,
                                visible=False),

                      FloatText(description='Nr. points per dimension',                                 # 12
                                value=10,
                                visible=False),

                      Dropdown(description='Consolidation Method',                                      # 13
                               values=['Average',
                                       'Geometric Mean'],
                               options=['Average',
                                        'Geometric Mean'],
                               value='Geometric Mean',
                               visible=False,
                               ),

                      Checkbox(description="Export Parameter Table (.xlsx)", value=False, visible=False)      # 14

                      ]

        column = VBox(children=new_column)
        column.header = new_column[0]
        column.vol_method = new_column[1]
        column.limit_vertices = new_column[2]

        column.nr_vertices = new_column[3]
        column.operating_point_method = new_column[4]
        column.enforce_bounds_operating_point = new_column[5]
        column.asymmetrical_bounds = new_column[6]

        column.bounds = new_column[7]
        column.lower_bound_par_space = new_column[7].children[0]
        column.upper_bound_par_space = new_column[7].children[1]

        column.ignore_unbounded = new_column[8]
        column.log_coordinates = new_column[9]
        column.delta = new_column[10]
        column.lamb = new_column[11]
        column.nr_points_sampling = new_column[12]
        column.consolidation_method = new_column[13]
        column.export_parameter_table = new_column[14]

        children.append(column)
        b.column_block.children = children
        b.remove_column.visible = True
        column.header.on_trait_change(self.update_column_block, 'value')
        column.vol_method.on_trait_change(self.update_column_block, 'value')
        column.limit_vertices.on_trait_change(self.update_column_block, 'value')
        column.operating_point_method.on_trait_change(self.update_column_block, 'value')
        fun = lambda name, value: self.generate_asymmetrical_bounds(name, value, column=column)
        column.asymmetrical_bounds.on_trait_change(fun, 'value')

    def generate_asymmetrical_bounds(self, name, asymmetrical, column=None):

        # This function should update the column.bounds box depending on the value of asymmetrical bounds.
        # If the bounds are asymmetrical, it should generate a list of dictionaries containing information of the bounds,
        # if the bounds are symmetrical, the column.bounds box should only contain a slider for lower and upper boudns.
        ind_variables = self.controller.ds.independent_variables

        if asymmetrical is True:
            asymm_list = []
            children = []
            for variable in ind_variables:
                asymm_dic = {}
                asymm_dic['name'] = variable
                asymm_dic['low_bound'] = FloatSlider(description='Lower Bound (Log)',
                                                     value=-3,
                                                     min=-9,
                                                     max=0,
                                                     step=0.001)

                asymm_dic['upper_bound'] = FloatSlider(description='Upper Bound (Log)',
                                                       value=3,
                                                       min=0,
                                                       max=9,
                                                       step=0.001,
                                                      )
                asymm_list.append(asymm_dic)
                children.append(HTML(value='<b>Bounds for ' + asymm_dic['name'] + '</b>'))
                children.append(asymm_dic['low_bound'])
                children.append(asymm_dic['upper_bound'])
            column.bounds.children = children
            column.bounds_list = asymm_list

        else:
            children = [
                            Slider(description='Lower Bound Par. Space (Log)',
                                   value=-3,
                                   min=-9,
                                   max=-1,
                                   step=1,
                                   ),

                            Slider(description='Upper Bound Par. Space (Log)',
                                   value=3,
                                   min=1,
                                   max=9,
                                   step=1,
                                   ),
                        ]

            column.bounds.children = children
            column.lower_bound_par_space = children[0]
            column.upper_bound_par_space = children[1]


    def add_neighbors_filter(self, b):
        # if len(b.column_block.children) == 0:
        #     return
        controller = self.controller
        children = [i for i in b.filters.children]
        new_filter = [Dropdown(description='Column #',
                               values=[str(i + 1) for i in range(1, 2 + len(b.column_block.children))],
                               options=[str(i + 1) for i in range(1, 2 + len(b.column_block.children))],
                               value=str(2 + len(b.column_block.children))),
                      Dropdown(description='Condition',
                               values=['==', '>', '<', '>=', '<=', '!='],
                               options=['==', '>', '<', '>=', '<=', '!=']),
                      Text(description='Value'),
                      Button(description='X')]
        new_filter[3].on_click(self.remove_neighbors_filter)
        current = HBox(children=new_filter)
        new_filter[3].current = current
        new_filter[3].filters = b.filters
        children.append(current)
        b.filters.children = children

    def remove_neighbors_filter(self, b):
        controller = self.controller
        if len(b.filters.children) == 0:
            return
        children = [i for i in b.filters.children]
        for i, child in enumerate(children):
            if b.current == child:
                children.pop(i)
        b.filters.children = children

    def save_table(self, b):
        controller = self.controller
        html_string = b.table_data
        controller.tables.add_table(html_string)
        controller.save_widget_data(b)

    def process_user_specified(self, l):
        cases = l.split(',')
        list_cases = []
        for case in cases:
            case_process = case.strip()
            if self.controller.ds(str(case_process)).is_valid() is True:
                list_cases.append(str(case_process))
        return list_cases

    def process_vertex_signature_neighbors(self, valid_phenotypes, b, case):

        controller = self.controller
        cases_list = []

        for c_string in valid_phenotypes:

            c = controller.ds(c_string)
            nr_shared_boundaries = len(case.shared_boundaries_indices(c, intersecting=b.intersecting.value)) if \
                case.shared_boundaries_indices(c, intersecting=b.intersecting.value) is not None else 0
            is_intersecting = case.share_boundaries_with(c, intersecting=b.intersecting.value)
            if case.shared_boundaries_is_valid(c) is False:
                continue
            if b.explicitly_share.value is True:
                if b.intersecting.value is True:
                    if nr_shared_boundaries != 0 and is_intersecting is True:
                        cases_list.append(c_string)
                else:
                    if nr_shared_boundaries != 0:
                        cases_list.append(c_string)

            else:  # they don't have to share a boundary
                if b.intersecting.value is True:
                    if is_intersecting is True:
                        cases_list.append(c_string)
                    else:
                        cases_list.append(c_string)

        return cases_list

    def calculate_all_mutation_rates_multiple(self):

        block = self.phenotypic_neighbors.all_mutation_rates
        numberOfMutations = block['nr_mutations'].value                                                 #0
        operating_point_method = block['operating_point_calculation_method'].value                      #1
        volume_method = block['volume_calculation_method'].value                                        #2
        enforce_bounds_operating_point = block['enforce_bounds_operating_point'].value                  #5
        enforce_bounds_volume = block['enforce_bounds_volume'].value                                    #6
        lower_bound = block['lower_bound_par_space'].value                                              #7
        upper_bound = block['upper_bound_par_space'].value                                              #8
        delta = block['delta_value'].value                                                              #9
        lamb = block['lambda_value'].value                                                              #10
        file_name = block['file_name'].value                                                            #13
        controller = self.controller
        par_values_volumes = {}
        par_values_k = {}
        par_values_k_vol = {}

        valid_cases = controller.ds.valid_cases()

        par_vals_dic = dict(('g' + str(ii), 1) for ii in valid_cases)
        par_vals_dic.update({'m': 1e-7, 's': 0.693, 'f': 1000})

        identity = controller.pidentity

        p_bounds = {}

        #####

        if enforce_bounds_operating_point is True:
            if block['asymmetrical_bounds_checkbox'].value is False:
                for key in controller.pvals.copy().keys():
                    p_bounds[key] = [10 ** lower_bound,  10 ** upper_bound]
            else:
                for el in block['asymmetrical_bounds_list']:
                    key = el['name']
                    lower_bound_i = el['low_bound'].value
                    upper_bound_i = el['upper_bound'].value
                    p_bounds[key] = [10 ** lower_bound_i, 10 ** upper_bound_i]
        else:
            p_bounds = None

        lb = {}
        ub = {}

        if enforce_bounds_volume is True:
            if block['asymmetrical_bounds_checkbox'].value is False:
                for key in controller.pvals.copy().keys():
                    lb[key] = 10 ** lower_bound
                    ub[key] = 10 ** upper_bound
            else:
                for el in block['asymmetrical_bounds_list']:
                    key = el['name']
                    lower_bound_i = el['low_bound'].value
                    upper_bound_i = el['upper_bound'].value
                    lb[key] = 10 ** lower_bound_i
                    ub[key] = 10 ** upper_bound_i
        else:
            for key in controller.pvals.copy().keys():
                lb[key] = 10 ** -20
                ub[key] = 10 ** 20

        #####

        par_values_normalization_factors = {}
        for from_case in valid_cases:
            key_from_normalization = str(from_case)
            from_case = controller.ds(from_case)
            if 'v'+from_case.case_number not in par_values_volumes:
                vol = from_case.volume(ignore_unbounded=False,
                                  log_coordinate=True,
                                  shared_boundaries=False,
                                  method=volume_method,
                                  # lowerBounds=lb,
                                  # upperBounds=ub,
                                  p_bounds=p_bounds,
                                  maxVertices=0,
                                  limitVertices=False)
                par_values_volumes['v'+from_case.case_number] = vol

            from_normalization = 0
            for to_case in valid_cases:
                to_case = controller.ds(to_case)
                if 'v' + to_case.case_number not in par_values_volumes:
                    vol = to_case.volume(ignore_unbounded=False,
                                           log_coordinate=True,
                                           shared_boundaries=False,
                                           method=volume_method,
                                           # lowerBounds=lb,
                                           # upperBounds=ub,
                                           p_bounds=p_bounds,
                                           maxVertices=0,
                                           limitVertices=False)
                    par_values_volumes['v' + to_case.case_number] = vol
                mutation_rate = from_case.mutation_rate_to_phenotype(to_case,
                                                                     identity,
                                                                     delta=delta,
                                                                     lamb=lamb,
                                                                     method=operating_point_method,
                                                                     p_bounds=p_bounds)
                key = from_case.mutation_key_to_case(to_case)
                par_values_k[key] = mutation_rate
                par_values_k_vol[key] = mutation_rate * par_values_volumes['v'+to_case.case_number]
                from_normalization += par_values_k_vol[key]
            par_values_normalization_factors['from_' + key_from_normalization] = from_normalization

        par_values_k_norm = {}
        for from_case in valid_cases:
            norm_factor = par_values_normalization_factors['from_' + str(from_case)]
            for to_case in valid_cases:
                key = mutation_key_to_case(from_case, to_case)
                par_values_k_norm[key] = par_values_k[key]/norm_factor

        par_vals_dic.update(par_values_k_norm)
        par_vals_dic.update(par_values_volumes)

        df = pd.DataFrame.from_dict(par_vals_dic, orient='index')
        df.columns = ['Value']
        df = df.rename_axis('Parameter', axis=1)
        df = df.sort_index()
        df.to_excel(file_name, index_label='Parameter')

    def get_mutation_rates_for_orthogonal_slice_ND(self, controller,
                                                   valid_cases_track,
                                                   constraints,
                                                   orthogonal_box,
                                                   p_bounds,
                                                   track_variable,
                                                   range_slice,
                                                   lamb, delta, identity,
                                                   lb, ub,
                                                   fractional_volumes_dic,
                                                   volume_method='Vertex Enumeration',
                                                   print_bol=False,
                                                   boundaries_for_track_variable=None,
                                                   tracks_method='Bounding Box',
                                                   width_method='Average of Methods',
                                                   include_fractional_volume=False,
                                                   include_lrs_volume=False,
                                                   ):

        mutation_rates_dic = {}

        # Let's calculate the bounding boxes for all phenotypes and the volumes for the alternative method.
        bounding_box_dic = {}
        volume_dic = {}
        lrs_volume = {}
        shape_dic = {}

        for case_nr in valid_cases_track:

            case = controller.ds(case_nr, constraints=constraints)
            bounding_box_dic[str(case_nr)] = case.bounding_box(p_bounds=p_bounds, log_out=True)

            if volume_method != 'Track Method':
                volume_dic[case_nr] = case.volume(ignore_unbounded=False,
                                                          log_coordinate=True,
                                                          shared_boundaries=False,
                                                          maxVertices=0,
                                                          limitVertices=False,
                                                          # lowerBounds=lb,
                                                          # upperBounds=ub,
                                                          p_bounds=p_bounds,
                                                          method=volume_method)
            else:
                vol, shape = get_track_volume_for_case_decomposing(controller, case_nr, track_variable,
                                                                                     constraints, range_slice,
                                                                                     boundaries_for_track_variable,
                                                                                     method=width_method,
                                                                                     p_bounds=p_bounds,
                                                                                     decomposing=False,
                                                                                     tracks_method=tracks_method)
                volume_dic[case_nr] = vol
                shape_dic[case_nr] = shape
            if include_lrs_volume is True:
                lrs_volume[case_nr] = case.volume(ignore_unbounded=False,
                                                  log_coordinate=True,
                                                  shared_boundaries=False,
                                                  maxVertices=0,
                                                  limitVertices=False,
                                                  # lowerBounds=lb,
                                                  # upperBounds=ub,
                                                  p_bounds=p_bounds,
                                                  method='Vertex Enumeration')

        # Main loop of the algorithm
        for from_case_key in valid_cases_track:

            from_case = controller.ds(from_case_key, constraints=constraints)
            from_case_OperatingPoint = from_case.valid_interior_parameter_set(p_bounds=p_bounds)
            vol_from = volume_dic[from_case_key]

            # ## Now we make some assessments
            if include_fractional_volume is True:

                    pvals_from_case = from_case_OperatingPoint.copy()
                    pvals_to_case = from_case_OperatingPoint.copy()
                    height_zero, height = is_height_zero(orthogonal_box, bounding_box_dic, from_case_key, from_case_key,
                                                         pvals_from_case, pvals_to_case)

                    tol_from_case = from_case.vertices_1D_slice(pvals_from_case,
                                                                track_variable,
                                                                range_slice=range_slice[track_variable],
                                                                log_out=True)
                    if len(tol_from_case) == 0:
                        continue

                    width_from = max(tol_from_case)[0] - min(tol_from_case)[0]
                    if width_from < 1.0e-10:
                        continue

                    fractional_volumes_dic['Box'].append(constraints)
                    fractional_volumes_dic['Phenotype'].append(from_case_key)
                    if from_case_key in shape_dic.keys():
                        fractional_volumes_dic['Formula'].append(shape_dic[from_case_key])
                    else:
                        fractional_volumes_dic['Formula'].append('NA')
                    if include_lrs_volume is True:
                        vol_from_lrs = lrs_volume[from_case_key]
                        fractional_volumes_dic['lrs Volume'].append(vol_from_lrs)
                        fractional_volumes_dic['Accurate'].append(round(vol_from_lrs, 6) == round(vol_from, 6))
                        fractional_volumes_dic['Absolute Error'].append(vol_from_lrs - vol_from)
                        fractional_volumes_dic['Relative Error'].append((vol_from_lrs - vol_from) / vol_from_lrs)
                    fractional_volumes_dic[volume_method + ' Volume'].append(vol_from)

                    width_from = max(tol_from_case)[0] - min(tol_from_case)[0]
                    centroid = min(tol_from_case)[0] + 0.5 * width_from
                    fractional_volumes_dic['Centroid Track Variable (log)'].append(centroid)

            # ## end assessments section

            for to_case_key in valid_cases_track:

                    to_case = controller.ds(to_case_key, constraints=constraints)
                    to_case_OperatingPoint = to_case.valid_interior_parameter_set(p_bounds=p_bounds)

                    key1 = from_case.mutation_key_to_case(to_case)
                    key2 = to_case.mutation_key_to_case(from_case)
                    if key1 in mutation_rates_dic.keys() or key2 in mutation_rates_dic.keys():
                        continue

                    pvals_from_case = from_case_OperatingPoint.copy()
                    pvals_to_case = to_case_OperatingPoint.copy()

                    mut_prob_from_to = 0
                    mut_prob_to_from = 0
                    height_zero, height = is_height_zero(orthogonal_box, bounding_box_dic, from_case_key, to_case_key,
                                                 pvals_from_case, pvals_to_case)
                    if height_zero is True:
                        continue

                    # now let's get the width for each phenotype.
                    # Here we need to calculate the tolerance and get the centroid value for the track parameter
                    tol_from_case = from_case.vertices_1D_slice(pvals_from_case,
                                                                track_variable,
                                                                range_slice=range_slice[track_variable],
                                                                log_out=True)
                    tol_to_case = to_case.vertices_1D_slice(pvals_to_case,
                                                            track_variable,
                                                            range_slice=range_slice[track_variable],
                                                            log_out=True)

                    if len(tol_from_case) == 0 or len(tol_to_case) == 0:
                        continue

                    width_from = max(tol_from_case)[0] - min(tol_from_case)[0]
                    width_to = max(tol_to_case)[0] - min(tol_to_case)[0]

                    if width_from < 1.0e-10 or width_to < 1.0e-10:
                        continue


                    pvals_from_case[track_variable] = 10**(min(tol_from_case)[0] + 0.5 * width_from)
                    pvals_to_case[track_variable] = 10**(min(tol_to_case)[0] + 0.5 * width_to)

                    if pvals_from_case[track_variable] == pvals_to_case[track_variable] and from_case_key != to_case_key:
                        continue

                    vol_to = volume_dic[to_case_key]

                    mut_prob_from_to += vol_to * DSPopDynamicsMutationRateForTransition(pvals_from_case._swigwrapper,
                                                                                        pvals_to_case._swigwrapper,
                                                                                        lamb,
                                                                                        delta,
                                                                                        identity._swigwrapper)

                    mut_prob_to_from += vol_from * DSPopDynamicsMutationRateForTransition(pvals_to_case._swigwrapper,
                                                                                          pvals_from_case._swigwrapper,
                                                                                          lamb,
                                                                                          delta,
                                                                                          identity._swigwrapper)
                    key = from_case.mutation_key_to_case(to_case, symbol='_')
                    mutation_rates_dic[key] = mut_prob_from_to
                    key = to_case.mutation_key_to_case(from_case, symbol='_')
                    mutation_rates_dic[key] = mut_prob_to_from


        return mutation_rates_dic

    def get_mutation_rates_for_orthogonal_slice_2D(self, controller,
                                                valid_cases_track,
                                                constraints, p_bounds,
                                                orthogonal_variable,
                                                track_variable,
                                                range_slice,
                                                lamb, delta, identity,
                                                lb, ub,
                                                volume_method='lrs'):

        mutation_rates_dic = {}
        par_copy = controller.pvals.copy()

        # Let's calculate the bounding boxes for all phenotypes
        bounding_box_dic = {}
        for case_nr in valid_cases_track:
            case = controller.ds(case_nr, constraints=constraints)
            bounding_box_dic[str(case_nr)] = case.bounding_box(p_bounds=p_bounds, log_out=True)

        # Main loop of the algorithm
        for from_case_key in valid_cases_track:
            from_case = controller.ds(from_case_key, constraints=constraints)
            from_case_OperatingPoint = from_case.valid_interior_parameter_set(p_bounds=p_bounds)

            for to_case_key in valid_cases_track:
                to_case = controller.ds(to_case_key, constraints=constraints)
                mut_prob_from_to = 0
                mut_prob_to_from = 0

                key1 = from_case.mutation_key_to_case(to_case)
                key2 = to_case.mutation_key_to_case(from_case)

                if key1 in mutation_rates_dic.keys() or key2 in mutation_rates_dic.keys():
                    continue

                to_case_OperatingPoint = to_case.valid_interior_parameter_set(p_bounds=p_bounds)

                count = 0
                width_from = 0
                width_to = 0
                height = 1
                orthogonal_bounds_dictionary = {}

                pvals_from_case = from_case_OperatingPoint.copy()
                pvals_to_case = to_case_OperatingPoint.copy()

                # We calculate the height for the orthogonal parameter
                height_i = min(bounding_box_dic[from_case_key][orthogonal_variable][1],
                               bounding_box_dic[to_case_key][orthogonal_variable][1]) - \
                           max(bounding_box_dic[from_case_key][orthogonal_variable][0],
                               bounding_box_dic[to_case_key][orthogonal_variable][0])

                # If height is negative or zero, then go to the next track parameter
                if height_i <= 0:
                    raise Exception("Height is zero")

                orthogonal_bounds = [max(bounding_box_dic[from_case_key][orthogonal_variable][0],
                                         bounding_box_dic[to_case_key][orthogonal_variable][0]),
                                     min(bounding_box_dic[from_case_key][orthogonal_variable][1],
                                         bounding_box_dic[to_case_key][orthogonal_variable][1]),
                                    ]

                orthogonal_coordinate = max(bounding_box_dic[from_case_key][orthogonal_variable][0],
                                            bounding_box_dic[to_case_key][orthogonal_variable][0]) + height_i * 0.5
                pvals_from_case[orthogonal_variable] = 10**orthogonal_coordinate

                pvals_to_case[orthogonal_variable] = 10**orthogonal_coordinate


                # now let's get the width for each phenotype.
                # Here we need to calculate the tolerance and get the centroid value for the track parameter

                tol_from_case = from_case.vertices_1D_slice(pvals_from_case,
                                                            track_variable,
                                                            range_slice=range_slice,
                                                            log_out=True)
                tol_to_case = to_case.vertices_1D_slice(pvals_to_case,
                                                        track_variable,
                                                        range_slice=range_slice,
                                                        log_out=True)

                if len(tol_from_case) == 0 or len(tol_to_case) == 0:
                    break

                # width_from_i = max(tol_from_case)[0] - min(tol_from_case)[0]
                # width_to_i = max(tol_to_case)[0] - min(tol_to_case)[0]

                #### testing modification on widths ###

                width_from_i = average_width(from_case,
                                             range_slice,
                                             pvals_from_case,
                                             track_variable,
                                             orthogonal_variable,
                                             orthogonal_bounds)

                width_to_i = average_width(to_case,
                                           range_slice,
                                           pvals_to_case,
                                           track_variable,
                                           orthogonal_variable,
                                           orthogonal_bounds)


                ### end testing

                width_from += width_from_i
                width_to += width_to_i
                height = height*height_i
                count += 1

                pvals_from_case[track_variable] = 10**(min(tol_from_case)[0] + 0.5*width_from_i)
                pvals_to_case[track_variable] = 10**(min(tol_to_case)[0] + 0.5*width_to_i)
                orthogonal_bounds_dictionary[orthogonal_variable] = orthogonal_bounds

                if count == 0:
                    continue

                A_width_from = width_from/count
                A_width_to = width_to/count
                vol_from = A_width_from * height
                vol_to = A_width_to * height

                ## Now we calculate real volumes for donor and recipient phenotypes
                #
                if len(orthogonal_bounds_dictionary.keys()) != 0 and volume_method != 'Track Method':   #and from_case.mutation_key_to_case(to_case) != 'k0315'
                    constraints_low = [str(key) + ' > ' + str(10**orthogonal_bounds_dictionary[key][0]) for key in
                                       orthogonal_bounds_dictionary.keys()]
                    constraints_up = [str(key) + ' < ' + str(10 ** orthogonal_bounds_dictionary[key][1]) for key in
                                       orthogonal_bounds_dictionary.keys()]
                    constraints = constraints_low + constraints_up

                    from_case_subset = controller.ds(from_case_key, constraints=constraints)
                    to_case_subset = controller.ds(to_case_key, constraints=constraints)

                    vol_from_real = from_case_subset.volume(ignore_unbounded=False,
                                                       log_coordinate=True,
                                                       shared_boundaries=False,
                                                       maxVertices=0,
                                                       limitVertices=False,
                                                       lowerBounds=lb,
                                                       upperBounds=ub,
                                                       method=volume_method)

                    vol_to_real = to_case_subset.volume(ignore_unbounded=False,
                                                   log_coordinate=True,
                                                   shared_boundaries=False,
                                                   maxVertices=0,
                                                   limitVertices=False,
                                                   lowerBounds=lb,
                                                   upperBounds=ub,
                                                   method=volume_method)


                # print("from_case, to_case: ", from_case_key, to_case_key)
                # print("pvals_from: ", pvals_from_case)
                # print("pvals_to: ", pvals_to_case)
                # print("volume from: ", vol_from_real, '==', vol_from)
                # print("vol to:", vol_to_real, '==', vol_to)

                mut_prob_from_to += vol_to * DSPopDynamicsMutationRateForTransition(pvals_from_case._swigwrapper,
                                                                                    pvals_to_case._swigwrapper,
                                                                                    lamb,
                                                                                    delta,
                                                                                    identity._swigwrapper)

                mut_prob_to_from += vol_from * DSPopDynamicsMutationRateForTransition(pvals_to_case._swigwrapper,
                                                                                      pvals_from_case._swigwrapper,
                                                                                      lamb,
                                                                                      delta,
                                                                                      identity._swigwrapper)
                key = from_case.mutation_key_to_case(to_case, symbol='_')
                mutation_rates_dic[key] = mut_prob_from_to
                key = to_case.mutation_key_to_case(from_case, symbol='_')
                mutation_rates_dic[key] = mut_prob_to_from

        return mutation_rates_dic

    def get_mutation_rates_over_track_variables_2D(self, controller,
                                                lb, ub, p_bounds,
                                                valid_cases,
                                                orthogonal_track_dictionary, lamb, delta, identity):

        variables = controller.pvals.copy().keys()
        mutation_rates_for_track_variables = {}
        range_slice = [lb, ub]

        for track_variable in variables:
            mutation_rates_for_orthogonal_variable = {}
            for orthogonal_variable in variables:
                if orthogonal_variable == track_variable:
                    continue

                # print("_________________________")
                # print("track_variable, orthogonal_variable: ", track_variable, orthogonal_variable)

                orthogonal_vertices = orthogonal_track_dictionary[orthogonal_variable]
                mutation_rates_for_orthogonal_tracks = []

                for i, orth_track_lb in enumerate(orthogonal_vertices):
                        if i + 1 == len(orthogonal_vertices):
                            break

                        orth_track_ub = orthogonal_vertices[i+1]
                        valid_cases_track = []
                        constraints = ''

                        for case_str in valid_cases:

                            constraints_low = [str(orthogonal_variable) + ' > ' + str(10**orth_track_lb)]
                            constraints_up = [str(orthogonal_variable) + ' < ' + str(10**orth_track_ub)]
                            constraints = constraints_low + constraints_up
                            case = controller.ds(case_str, constraints=constraints)
                            if case.is_valid(p_bounds=p_bounds) is True:
                                valid_cases_track.append(case.case_number)

                        # Now se send the subset of valid cases "valid_cases_track" to a routine to calculate a dictionary
                        # of mutation rates

                        if len(valid_cases_track) != 0:
                            # print("orthogonal track: ", i+1)
                            mut_rates_dic = self.get_mutation_rates_for_orthogonal_slice_2D(controller,
                                                                                    valid_cases_track,
                                                                                    constraints,
                                                                                    p_bounds,
                                                                                    orthogonal_variable,
                                                                                    track_variable,
                                                                                    range_slice,
                                                                                    lamb, delta, identity, lb, ub
                                                                                    )
                        else:
                            mut_rates_dic = {}
                        mutation_rates_for_orthogonal_tracks.append(mut_rates_dic)

                        # print("track variable, orthogonal variable, lb, ub, valid cases:", track_variable,
                        #       orthogonal_variable, orth_track_lb, orth_track_ub, valid_cases_track)
                mutation_rates_for_orthogonal_variable[orthogonal_variable] = mutation_rates_for_orthogonal_tracks
            mutation_rates_for_track_variables[track_variable] = mutation_rates_for_orthogonal_variable

        return mutation_rates_for_track_variables

    def get_mutation_rates_over_track_variables_ND(self, controller,
                                                   lb, ub, p_bounds,
                                                   valid_cases,
                                                   orthogonal_track_dictionary, lamb, delta, identity, volume_method,
                                                   tracks_method='Bounding Box',
                                                   width_method='Average of Methods',
                                                   include_fractional_volume=False,
                                                   include_lrs_volume=False):

        variables = controller.pvals.copy().keys()
        mutation_rates_for_track_variables = {}
        # range_slice = [lb, ub]
        range_slice = p_bounds
        orthogonal_track_pair_dictionary = process_orthogonal_track_dictionary(orthogonal_track_dictionary)
        fractional_volumes_for_track_variable_dic = {}

        for track_variable in variables:

            orthogonal_boxes = get_orthogonal_box(track_variable, orthogonal_track_pair_dictionary)

            # We now expand the list of orthogonal boxes to include any potential boxes
            orthogonal_boxes_expanded = get_additional_boxes(track_variable,
                                                             orthogonal_boxes,
                                                             valid_cases,
                                                             p_bounds,
                                                             controller, lb, ub, method=tracks_method)

            mutation_rates_for_track_variable = []

            if include_lrs_volume is True:
                fractional_volumes_dic = {'Box': [], 'Phenotype': [], volume_method + ' Volume': [], 'lrs Volume': [],
                                          'Centroid Track Variable (log)': [], 'Accurate': [], 'Absolute Error': [],
                                          'Relative Error': [], 'Formula': []}
            else:
                fractional_volumes_dic = {'Box': [], 'Phenotype': [], volume_method + ' Volume': [],
                                          'Centroid Track Variable (log)': [], 'Formula': []}

            for orthogonal_box in orthogonal_boxes_expanded:

                        valid_cases_track = []
                        constraints, vol_box = get_constraint_list(orthogonal_box)

                        for case_str in valid_cases:

                                case = controller.ds(case_str, constraints=constraints)
                                if case.is_valid(p_bounds=p_bounds) is True:
                                        # now we make sure that the phenotype is not a dot or a line.
                                        skip = False
                                        bb_case = case.bounding_box(log_out=True)

                                        for key in bb_case.keys():
                                                if len(bb_case[key]) == 1 or len(bb_case[key]) == 0:
                                                    skip = True
                                                    break
                                                if max(bb_case[key]) - min(bb_case[key]) < 1e-10:
                                                    skip = True
                                                    break

                                        if skip is False:
                                            valid_cases_track.append(case.case_number)

                        boundaries_for_track_variable = get_boundaries_for_track_variable(controller,
                                                                                          valid_cases_track,
                                                                                          lb, ub, p_bounds,
                                                                                          track_variable,
                                                                                          constraints=constraints,
                                                                                          method=tracks_method
                                                                                          )

                        # Now se send the subset of valid cases "valid_cases_track" to a routine to calculate a dictionary
                        # of mutation rates
                        if len(valid_cases_track) != 0:
                            mut_rates_dic = self.get_mutation_rates_for_orthogonal_slice_ND(controller,
                                                                                            valid_cases_track,
                                                                                            constraints,
                                                                                            orthogonal_box,
                                                                                            p_bounds,
                                                                                            track_variable,
                                                                                            range_slice,
                                                                                            lamb, delta, identity,
                                                                                            lb, ub,
                                                                                            fractional_volumes_dic,
                                                                                            volume_method=volume_method,
                                                                                            boundaries_for_track_variable=boundaries_for_track_variable,
                                                                                            tracks_method=tracks_method,
                                                                                            width_method=width_method,
                                                                                            include_fractional_volume=include_fractional_volume,
                                                                                            include_lrs_volume=include_lrs_volume,
                                                                                            )
                        else:
                            mut_rates_dic = {}
                        mutation_rates_for_track_variable.append(mut_rates_dic)
            mutation_rates_for_track_variables[track_variable] = mutation_rates_for_track_variable
            fractional_volumes_for_track_variable_dic[track_variable] = fractional_volumes_dic

        return mutation_rates_for_track_variables, fractional_volumes_for_track_variable_dic

    def get_normalized_rates_ND(self, controller, mut_rates_over_track_variable, valid_cases, symbol='_'):

        # Method 1. Normalization factors are calculated over all track variables.

        # normalized_mutation_rates_dic = {}
        # normalization_factors_from = {}
        # variables = mut_rates_over_track_variable.keys()

        # # We first create the normalization factors
        # for from_case_key in valid_cases:
        #
        #     factor = 0
        #     key = 'from_' + str(from_case_key)
        #
        #     for track_variable in variables:
        #         for orthogonal_box in mut_rates_over_track_variable[track_variable]:
        #
        #             for transition in orthogonal_box.keys():
        #                     if int(from_case_key) == int(transition.split('_')[0].split('k')[1]):
        #                         factor += orthogonal_box[transition]
        #
        #     normalization_factors_from[key] = factor
        #
        # # to now create the normalized rates
        # for from_case_key in valid_cases:
        #     key_from = 'from_' + str(from_case_key)
        #
        #     for to_case_key in valid_cases:
        #
        #         key = mutation_key_to_case(from_case_key, to_case_key, symbol=None)
        #         key_ = mutation_key_to_case(from_case_key, to_case_key, symbol="_")
        #         mutation_rate = 0
        #
        #         for track_variable in variables:
        #
        #             for orthogonal_box in mut_rates_over_track_variable[track_variable]:
        #                 if key_ in orthogonal_box.keys():
        #                     mutation_rate += orthogonal_box[key_]
        #
        #         normalized_mutation_rates_dic[key] = mutation_rate / normalization_factors_from[key_from]

        # Method 2. Normalization factors are created for each track variable

        normalized_mutation_rates_dic = {}
        normalization_factors_from_over_track_variables = {}
        variables = mut_rates_over_track_variable.keys()
        normalized_mutation_rates_track_variable_dic = {}

        # We first create the normalization factors

        for track_variable in variables:

            normalization_factors_from = {}
            for from_case_key in valid_cases:
                factor = 0
                key = 'from_' + str(from_case_key)
                for orthogonal_box in mut_rates_over_track_variable[track_variable]:
                    for transition in orthogonal_box.keys():
                            if int(from_case_key) == int(transition.split('_')[0].split('k')[1]):
                                factor += orthogonal_box[transition]

                    normalization_factors_from[key] = factor
            normalization_factors_from_over_track_variables[track_variable] = normalization_factors_from

        # to now create the normalized rates
        for from_case_key in valid_cases:
            key_from = 'from_' + str(from_case_key)

            for to_case_key in valid_cases:

                key = mutation_key_to_case(from_case_key, to_case_key, symbol=None)
                key_ = mutation_key_to_case(from_case_key, to_case_key, symbol="_")
                mutation_rate = 0

                for track_variable in variables:

                    for orthogonal_box in mut_rates_over_track_variable[track_variable]:
                        if key_ in orthogonal_box.keys():
                            mutation_rate += orthogonal_box[key_] / normalization_factors_from_over_track_variables[track_variable][key_from]

                normalized_mutation_rates_dic[key] = mutation_rate

        # now calculate normalized mutation rates for each track variable
        for track_variable in variables:
            mut_rates_dic = {}

            for from_case_key in valid_cases:
                key_from = 'from_' + str(from_case_key)

                for to_case_key in valid_cases:

                    key = mutation_key_to_case(from_case_key, to_case_key, symbol=None)
                    key_ = mutation_key_to_case(from_case_key, to_case_key, symbol="_")
                    mutation_rate = 0

                    for orthogonal_box in mut_rates_over_track_variable[track_variable]:
                        if key_ in orthogonal_box.keys():
                            mutation_rate += orthogonal_box[key_] / normalization_factors_from_over_track_variables[track_variable][key_from]
                    mut_rates_dic[key] = mutation_rate
            normalized_mutation_rates_track_variable_dic[track_variable] = mut_rates_dic

        return normalized_mutation_rates_dic, normalized_mutation_rates_track_variable_dic

    def get_normalized_rates_2D(self, controller, mut_rates_over_track_variable, valid_cases, symbol='_'):

        normalized_mutation_rates_dic = {}
        normalization_factors_from = {}

        variables = mut_rates_over_track_variable.keys()

        # We first create the normalization factors
        for from_case_key in valid_cases:

            factor = 0
            key = 'from_' + str(from_case_key)

            for track_variable in variables:
                for orthogonal_variable in variables:
                    if orthogonal_variable == track_variable:
                        continue
                    for mutations_track in mut_rates_over_track_variable[track_variable][orthogonal_variable]:
                        for transition in mutations_track.keys():
                            if int(from_case_key) == int(transition.split('_')[0].split('k')[1]):
                                factor += mutations_track[transition]
            normalization_factors_from[key] = factor

        # to now create the normalized rates
        for from_case_key in valid_cases:
            key_from = 'from_' + str(from_case_key)

            for to_case_key in valid_cases:

                key = mutation_key_to_case(from_case_key, to_case_key, symbol=None)
                key_ = mutation_key_to_case(from_case_key, to_case_key, symbol="_")
                mutation_rate = 0

                for track_variable in variables:
                    for orthogonal_variable in variables:
                        if orthogonal_variable == track_variable:
                            continue
                        for mutations_track in mut_rates_over_track_variable[track_variable][orthogonal_variable]:
                            if key_ in mutations_track.keys():
                                mutation_rate += mutations_track[key_]
                normalized_mutation_rates_dic[key] = mutation_rate/normalization_factors_from[key_from]

        return normalized_mutation_rates_dic

    def calculate_all_mutation_rates_single_ND(self):

        block = self.phenotypic_neighbors.all_mutation_rates
        volume_method = block['volume_calculation_method'].value                                            #2
        tracks_method = block['tracks_calculation_method'].value                                            #3
        width_method = block['method_for_width_determination'].value                                        #4
        enforce_bounds_operating_point = block['enforce_bounds_operating_point'].value                      #5
        enforce_bounds_volume = block['enforce_bounds_volume'].value                                        #6
        lower_bound = block['lower_bound_par_space'].value                                                  #7
        upper_bound = block['upper_bound_par_space'].value                                                  #8
        delta = block['delta_value'].value                                                                  #9
        lamb = block['lambda_value'].value                                                                  #10
        include_fractional_volume = block['include_fractional_volume_data'].value                           #11
        include_lrs_volume = block['include_real_volumes'].value                                            #12
        file_name = block['file_name'].value                                                                #13
        controller = self.controller
        valid_cases = controller.ds.valid_cases()
        identity = controller.pidentity

        par_copy = controller.pvals.copy()
        p_bounds = {}

        if identity is None:
            raise ValueError("Please set the identity of each parameter in the Edit Parameters tab")

        par_vals_dic = dict(('g' + str(ii), 1) for ii in valid_cases)
        par_vals_dic.update({'m': 1e-7, 's': 0.693, 'f': 1000})

        if enforce_bounds_operating_point is True:
            # range_slice = [10 ** lower_bound,  10 ** upper_bound]

            if block['asymmetrical_bounds_checkbox'].value is False:
                for key in par_copy.keys():
                    p_bounds[key] = [10 ** lower_bound,  10 ** upper_bound]
            else:
                for el in block['asymmetrical_bounds_list']:
                    key = el['name']
                    lower_bound_i = el['low_bound'].value
                    upper_bound_i = el['upper_bound'].value
                    p_bounds[key] = [10 ** lower_bound_i, 10 ** upper_bound_i]
        else:
            p_bounds = None
            # range_slice = None

        lb = {}
        ub = {}

        if enforce_bounds_volume is True:
            if block['asymmetrical_bounds_checkbox'].value is False:
                for key in par_copy.keys():
                    lb[key] = 10 ** lower_bound
                    ub[key] = 10 ** upper_bound
            else:
                for el in block['asymmetrical_bounds_list']:
                    key = el['name']
                    lower_bound_i = el['low_bound'].value
                    upper_bound_i = el['upper_bound'].value
                    lb[key] = 10 ** lower_bound_i
                    ub[key] = 10 ** upper_bound_i
        else:
            for key in par_copy.keys():
                lb[key] = 10 ** -20
                ub[key] = 10 ** 20

        orthogonal_track_dictionary = get_track_dictionary(controller, valid_cases, lb, ub, p_bounds,
                                                           method=tracks_method)

        mut_rates_over_track_variable, \
        fractional_volumes_for_track_variable_dic = self.get_mutation_rates_over_track_variables_ND(controller,
                                                                                        lb,
                                                                                        ub,
                                                                                        p_bounds,
                                                                                        valid_cases,
                                                                                        orthogonal_track_dictionary,
                                                                                        lamb, delta, identity,
                                                                                        volume_method,
                                                                                        tracks_method=tracks_method,
                                                                                        width_method=width_method,
                                                                                        include_fractional_volume=include_fractional_volume,
                                                                                        include_lrs_volume=include_lrs_volume)

        normalized_mutation_rates_dic, \
        normalized_mutation_rates_track_variable_dic = self.get_normalized_rates_ND(controller,
                                                                                     mut_rates_over_track_variable,
                                                                                     valid_cases,
                                                                                     symbol='_')
        # # export to Excel
        with pd.ExcelWriter(file_name) as writer:

            updated_dic = par_vals_dic.copy()
            updated_dic.update(normalized_mutation_rates_dic)
            df = pd.DataFrame.from_dict(updated_dic, orient='index')
            df.columns = ['Value']
            df = df.rename_axis('Parameter', axis=1)
            df = df.sort_index()
            df.to_excel(writer, sheet_name='All_Track_Variables') #, index_label='Parameter'

            for track_variable in normalized_mutation_rates_track_variable_dic.keys():
                updated_dic = par_vals_dic.copy()
                updated_dic.update(normalized_mutation_rates_track_variable_dic[track_variable])
                df = pd.DataFrame.from_dict(updated_dic, orient='index')
                df.columns = ['Value']
                df = df.rename_axis('Parameter', axis=1)
                df = df.sort_index()
                df.to_excel(writer, sheet_name='Track_Variable_ ' + track_variable, index_label='Parameter')

        # Now export to Excel fractional volumes
        if include_fractional_volume is True:

            if include_lrs_volume is True:
                order = ['Box', 'Phenotype', volume_method + ' Volume', 'Formula', 'lrs Volume', 'Centroid Track Variable (log)',
                         'Accurate', 'Absolute Error', 'Relative Error']
            else:
                order = ['Box', 'Phenotype', volume_method + ' Volume', 'Formula', 'Centroid Track Variable (log)']

            with pd.ExcelWriter(self.controller.name + '_Fractional_Volumes' +'.xlsx') as writer:
                for track_variable in fractional_volumes_for_track_variable_dic.keys():
                        df = pd.DataFrame(data=fractional_volumes_for_track_variable_dic[track_variable])
                        df = df[order]
                        df.to_excel(writer, sheet_name='Track_Variable_'+track_variable, index=False)

    def get_neighbors(self, b):
        # The target output is self.phenotypic_neighbors.children
        # Calculate Phenotypic neighbors

        if str(b.phenotype.value) == '' and b.neighbors_method.value != 'All Mutation Rates':
            return

        if b.description == 'Modify Table':
            # self.update_phenotypic_neighbors(bol_intersecting=b.intersecting.value, bol_share=b.explicitly_share.value,
            #                                  phen_value=b.phenotype.value, extra_columns=self.extra_columns,
            #                                  method=b.neighbors_method.value)
            self.update_phenotypic_neighbors()
            return

        if b.neighbors_method.value == 'All Mutation Rates':

            block = self.phenotypic_neighbors.all_mutation_rates
            if block['nr_mutations'].value == 'Multiple':
                self.calculate_all_mutation_rates_multiple()
            else:
                self.calculate_all_mutation_rates_single_ND()

            return

        controller = self.controller
        headers = []

        # 1. Get case number and valid phenotypes
        case_s = str(b.phenotype.value)
        case = controller.ds(case_s)

        valid_phenotypes = []

        if b.neighbors_method.value == 'Signature':
            valid_phenotypes = case.neighbors_signature(controller.ds)
            valid_phenotypes = self.process_vertex_signature_neighbors(valid_phenotypes, b, case)
        elif b.neighbors_method.value == 'Vertex Enumeration':
            valid_phenotypes = self.process_vertex_signature_neighbors(controller.ds.valid_cases(), b, case)
        elif b.neighbors_method.value == 'User Defined':
            valid_phenotypes = self.process_user_specified(b.user_cases.value)

        if len(valid_phenotypes) == 0:
            return

        s = '<div><table>\n<caption><b>Neighbors for Case {} in the System Design Space. </b></caption>\n'.format(case.case_number)
        s += '<tr align=center><td style="padding:0 15px 0 15px;"><b>{0}</b>' \
             '</td><td style="padding:0 15px 0 15px;"><b>{1}</b></td>'.format('  Case Number  ',
                                                                              '  Case Signature  ')
        headers += ['Case Number', 'Case Signature']

        for column in b.extra_columns.children:

            if column.header.value == 'Number of Vertices':
                header = 'Number <br> of Vertices'
            elif column.header.value == 'Dimension of Shared Space':
                header = 'Dimension <br> of Shared Space'
            elif column.header.value == 'Number of Shared Boundaries':
                header = 'Number <br> of Shared Boundaries'
            elif column.header.value == 'Volume of Shared Space':
                header = 'Volume <br> Shared Space'
                if column.vol_method.value == 'Tolerances':
                    header += ' (Tolerances)'
                else:
                    header += ' (lrs)'
            elif column.header.value == 'Distance':
                header = 'Distance'
            elif column.header.value == 'Intersect':
                header = 'Intersect'
            elif column.header.value == 'Mutation Rate':
                header = 'Mutation Rate'
            headers += [header]
            s += '<td><b>' + header + '</b></td>'
            column.intersecting = b.intersecting
        s += '</tr>'

        b.headers = headers
        d = dict((i, []) for i in b.headers)

        #2. Loop over valid cases and determine if cases intersect and share boundaries or not.
        for c_string in valid_phenotypes:

            if c_string == case_s:
                continue

            c = controller.ds(c_string)

            values = [c.case_number, c.signature]
            values += [self.value_for_extra_column(case, c, column) for column in
                       b.extra_columns.children]
            if b.export_table.value is True:
                for i in range(len(b.headers)):
                    d[b.headers[i]].append(values[i])
            if self.show_case(values, b.extra_columns.children, b.filters) is False:
                continue
            s += create_string_for_table(c)
            for column in b.extra_columns.children:
                s += self.html_for_extra_column(case, c, column)
            s += '</tr>\n'


        s += '</table><caption>'
        s += '<br>The number of neighboring cases listed in this table is: ' + str(s.count("</tr>\n"))
        s += '</caption></div>'
        b.table_data = d

        #3. Assign HTML table
        table = HTML(value=s)
        b.description = 'Modify Table'
        save_table = Button(description='Save Table')
        save_table.table_data = s
        save_table.on_click(self.save_table)

        if old_ipython is True:
            table.set_css('height', '300px')
        else:
            table.height = '300px'
            table.overflow_x = 'auto'
            table.overflow_y = 'auto'

        #4. Generate Excel file
        if b.export_table.value is True:
            df = pd.DataFrame(data=b.table_data)
            df.to_excel(controller.name + '_neighbors' + '.xlsx', index=False, columns=headers)

        self.phenotypic_neighbors.children = [b, save_table, table]

    def html_for_extra_column(self, case1, case2, column):

        controller = self.controller
        s = '<td style="padding:0 15px 0 15px;">'

        if column.header.value == 'Number of Vertices':
            lb = controller.pvals.copy()
            up = controller.pvals.copy()
            maxVertices = int(column.nr_vertices.value)
            for key in lb.keys():
                lb[key] = 10**column.lower_bound_par_space.value          #1E-6
                up[key] = 10**column.upper_bound_par_space.value          #1E6
            s += str(case1.shared_boundaries_number_of_vertices(case2, lb, up, maxVertices, column.limit_vertices.value))

        elif column.header.value == 'Dimension of Shared Space':
            case_int = dspace.CaseIntersection([case1, case2])
            lb = controller.pvals.copy()
            ub = controller.pvals.copy()
            maxVertices = int(column.nr_vertices.value)

            for key in lb.keys():
                lb[key] = 10**column.lower_bound_par_space.value          #1E-6
                ub[key] = 10**column.upper_bound_par_space.value          #1E6
            s += str(case_int.dimension(lb, ub, maxVertices))
        elif column.header.value == 'Number of Shared Boundaries':
            nr_shared_boundaries = len(case1.shared_boundaries_indices(case2, intersecting=column.intersecting.value))
            s += str(nr_shared_boundaries)
        elif column.header.value == 'Volume of Shared Space':
            case_int = dspace.CaseIntersection([case1, case2])
            lb = controller.pvals.copy()
            ub = controller.pvals.copy()
            maxVertices = int(column.nr_vertices.value)

            for key in lb.keys():
                lb[key] = 10**column.lower_bound_par_space.value          #1E-6
                ub[key] = 10**column.upper_bound_par_space.value          #1E6

            if column.vol_method.value == 'Vertex Enumeration':
                vol, _ = case_int.volume_lrs(lb, ub, maxVertices, column.limit_vertices.value)
            else:
                ignore_unbounded = column.ignore_unbounded.value
                log_coordinate = column.log_coordinates.value
                vol = case_int.volume(ignore_unbounded=ignore_unbounded,
                                      log_coordinate=log_coordinate,
                                      shared_boundaries=True)

            s += '%.2E' % Decimal(str(round_sig(vol)))

        elif column.header.value == 'Distance':
            dist = case1.calculate_distance_to_case(case2)
            s += '%.2E' % Decimal(str(round_sig(dist)))

        elif column.header.value == 'Intersect':
            case_int = dspace.CaseIntersection([case1, case2])
            inters = case_int.is_valid()
            s += '+' if inters is True else '-'
        elif column.header.value == 'Mutation Rate':

            identity = controller.pidentity
            delta = column.delta.value
            lamb = column.lamb.value
            nr_points = int(column.nr_points_sampling.value)
            consolidation_method = column.consolidation_method.value
            operating_point_method = column.operating_point_method.value
            enforce_bounds_operating_point = column.enforce_bounds_operating_point.value
            upper_bound = column.upper_bound_par_space.value
            lower_bound = column.lower_bound_par_space.value
            asymmetrical_bound = column.asymmetrical_bounds.value
            if asymmetrical_bound is True:
                dic_list = column.bounds_list
            p_bounds = {}
            for key in controller.pvals.copy().keys():
                if enforce_bounds_operating_point is True:
                    if asymmetrical_bound is False:
                        p_bounds[key] = [10 ** lower_bound, 10 ** upper_bound]
                    else:
                        for var_dic in dic_list:
                            p_bounds[var_dic['name']] = [10 ** var_dic['low_bound'].value,
                                                         10 ** var_dic['upper_bound'].value]
                else:
                    p_bounds = None

            mutation = case1.mutation_rate_to_phenotype(case2,
                                                        identity,
                                                        delta=delta,
                                                        lamb=lamb,
                                                        method=operating_point_method,
                                                        p_bounds=p_bounds,
                                                        nr_points=nr_points,
                                                        average_method=consolidation_method
                                                        )

            s += '%.2E' % Decimal(str(round_sig(mutation)))
        else:
            s += '-'
        s += '</td>'
        return s

    def value_for_extra_column(self, case1, case2, column):
        controller = self.controller
        if column.header.value == 'Number of Vertices':
            lb = controller.pvals.copy()
            ub = controller.pvals.copy()
            maxVertices = int(column.nr_vertices.value)

            for key in lb.keys():
                lb[key] = 10**column.lower_bound_par_space.value          #1E-6
                ub[key] = 10**column.upper_bound_par_space.value          #1E6
            value = case1.shared_boundaries_number_of_vertices(case2, lb, ub, maxVertices, column.limit_vertices.value)

        elif column.header.value == 'Dimension of Shared Space':
            case_int = dspace.CaseIntersection([case1, case2])
            lb = controller.pvals.copy()
            ub = controller.pvals.copy()
            maxVertices = int(column.nr_vertices.value)

            for key in lb.keys():
                lb[key] = 10**column.lower_bound_par_space.value              #1E-6
                ub[key] = 10**column.upper_bound_par_space.value              #1E6
            value = case_int.dimension(lb, ub, maxVertices)

        elif column.header.value == 'Number of Shared Boundaries':
            nr_shared_boundaries = len(case1.shared_boundaries_indices(case2, intersecting=column.intersecting.value))
            value = nr_shared_boundaries

        elif column.header.value == 'Volume of Shared Space':
            case_int = dspace.CaseIntersection([case1, case2])
            lb = controller.pvals.copy()
            ub = controller.pvals.copy()
            maxVertices = int(column.nr_vertices.value)

            for key in lb.keys():
                lb[key] = 10**column.lower_bound_par_space.value          #1E-6
                ub[key] = 10**column.upper_bound_par_space.value          #1E6

            if column.vol_method.value == 'Vertex Enumeration':
                value, _ = case_int.volume_lrs(lb, ub, maxVertices, column.limit_vertices.value)
            else:
                ignore_unbounded = column.ignore_unbounded.value
                log_coordinate = column.log_coordinates.value
                value = case_int.volume(ignore_unbounded=ignore_unbounded,
                                        log_coordinate=log_coordinate,
                                        shared_boundaries=True)

        elif column.header.value == 'Distance':
            value = case1.calculate_distance_to_case(case2)

        elif column.header.value == 'Intersect':
            inters = dspace.CaseIntersection([case1, case2])
            value = '+' if inters.is_valid() is True else '-'

        elif column.header.value == 'Mutation Rate':
            identity = controller.pidentity
            delta = column.delta.value
            lamb = column.lamb.value
            nr_points = int(column.nr_points_sampling.value)
            consolidation_method = column.consolidation_method.value
            operating_point_method = column.operating_point_method.value
            enforce_bounds_operating_point = column.enforce_bounds_operating_point.value
            upper_bound = column.upper_bound_par_space.value
            lower_bound = column.lower_bound_par_space.value
            asymmetrical_bound = column.asymmetrical_bounds.value
            if asymmetrical_bound is True:
                dic_list = column.bounds_list

            p_bounds = {}
            for key in controller.pvals.copy().keys():
                if enforce_bounds_operating_point is True:
                    if asymmetrical_bound is False:
                        p_bounds[key] = [10 ** lower_bound, 10 ** upper_bound]
                    else:
                        for var_dic in dic_list:
                            p_bounds[var_dic['name']] = [10 ** var_dic['low_bound'].value,
                                                         10 ** var_dic['upper_bound'].value]
                else:
                    p_bounds = None

            # mutation = case1.mutation_rate_to_phenotype(case2,
            #                                             identity,
            #                                             delta=delta,
            #                                             lamb=lamb,
            #                                             method=operating_point_method,
            #                                             p_bounds=p_bounds,
            #                                             nr_points=nr_points,
            #                                             average_method=consolidation_method
            #                                             )
            mutation = 1
            value = mutation

        else:
            value = 0

        return value

    def show_case(self, values, columns, filters):
        showCase = True
        for filter_box in filters.children:
            a_filter = filter_box.children
            column_number = int(a_filter[0].value)
            if column_number == 1:
                rhs = int(a_filter[2].value)
                lhs = int(values[column_number - 1])
            elif column_number == 2:
                rhs = str(a_filter[2].value).replace(" ", "")
                lhs = str(values[column_number - 1]).replace(" ", "")
                if rhs.find("*") != -1:
                    siglist = self.generate_multiple_signatures(rhs)
                    showCase = self.show_case_multiple_signatures(siglist, a_filter, lhs)
                    break
            else:
                column_type = columns[column_number - 3].header.value
                if column_type == 'Validity' or (column_type == '# eigenvalues w/ positive real part' and
                                                 columns[column_number - 3].root_method.value == 'Numpy Functions' and
                                                 columns[column_number - 3].show_complex_conjugates.value is True)\
                        or (column_type == 'Is a Blowup') or (column_type == 'Is False Blowing'):
                    rhs = str(a_filter[2].value)
                    lhs = str(values[column_number - 1])
                else:
                    if values[column_number - 1] == '-' or values[column_number - 1] == '*':
                        showCase = False
                        break
                    rhs = float(a_filter[2].value)
                    lhs = float((str(values[column_number - 1])).replace("*",
                                                                         ""))
            if a_filter[1].value == '==' and lhs != rhs:
                showCase = False
                break
            if a_filter[1].value == '>=' and lhs < rhs:
                showCase = False
                break
            if a_filter[1].value == '<=' and lhs > rhs:
                showCase = False
                break
            if a_filter[1].value == '>' and lhs <= rhs:
                showCase = False
                break
            if a_filter[1].value == '<' and lhs >= rhs:
                showCase = False
                break
            if a_filter[1].value == '!=' and lhs == rhs:
                showCase = False
                break
        return showCase

    def generate_multiple_signatures(self, rhs):

        signature = rhs
        controller = self.controller.ds
        i = 0
        siglist = []
        wild_cards = []
        while i < len(signature):
            if signature[i] == '(':
                start = i+1
                while signature[i] != ')':
                    i += 1
                siglist.append(int(signature[start:i]))
            elif signature[i] == '*':
                num_wild = int(controller.signature[i])
                for j in range(num_wild):
                    new_sig = signature.replace('*', str(j+1), 1)
                    wild_cards += self.generate_multiple_signatures(new_sig)
                return wild_cards
            else:
                siglist.append(int(signature[i]))
            i += 1

        return [siglist]

    def show_case_multiple_signatures(self, siglist, a_filter, lhs):
        showCase = True
        for signature_int in siglist:
            rhs = ''.join(str(x) for x in signature_int).replace(" ", "")
            if a_filter[1].value == '==' and lhs != rhs:
                showCase = False
                break
            if a_filter[1].value == '!=' and lhs == rhs:
                showCase = False
                break
        return showCase

    def update_column_block(self, name, value):

        for column in self.extra_columns.children:

            column.delta.visible = False
            column.lamb.visible = False
            column.operating_point_method.visible = False
            column.enforce_bounds_operating_point.visible = False
            column.asymmetrical_bounds.visible = False
            column.export_parameter_table.visible = False
            column.nr_points_sampling.visible = False
            column.consolidation_method.visible = False

            if column.header.value == 'Number of Vertices':
                column.limit_vertices.visible = True
                if column.limit_vertices.value is True:
                    column.nr_vertices.visible = True
                else:
                    column.nr_vertices.visible = False
                column.vol_method.visible = False
                column.upper_bound_par_space.visible = True
                column.lower_bound_par_space.visible = True

            elif column.header.value == 'Volume of Shared Space':
                column.vol_method.visible = True
                if column.vol_method.value == 'Vertex Enumeration':
                    column.limit_vertices.visible = True
                    column.upper_bound_par_space.visible = True
                    column.lower_bound_par_space.visible = True
                    if column.limit_vertices.value is True:
                        column.nr_vertices.visible = True
                    else:
                        column.nr_vertices.visible = False
                else:
                    column.limit_vertices.visible = False
                    column.nr_vertices.visible = False
                    column.upper_bound_par_space.visible = False
                    column.lower_bound_par_space.visible = False
            elif column.header.value == 'Dimension of Shared Space':
                column.limit_vertices.visible = True
                if column.limit_vertices.value is True:
                    column.nr_vertices.visible = True
                else:
                    column.nr_vertices.visible = False
                column.vol_method.visible = False
                column.upper_bound_par_space.visible = True
                column.lower_bound_par_space.visible = True
            elif column.header.value == 'Number of Shared Boundaries':
                column.limit_vertices.visible = False
                column.nr_vertices.visible = False
                column.vol_method.visible = False
                column.upper_bound_par_space.visible = False
                column.lower_bound_par_space.visible = False

            elif column.header.value == 'Mutation Rate':
                column.delta.visible = True
                column.lamb.visible = True
                column.vol_method.visible = True
                column.operating_point_method.visible = True
                column.enforce_bounds_operating_point.visible = True
                column.asymmetrical_bounds.visible = True
                column.upper_bound_par_space.visible = True
                column.lower_bound_par_space.visible = True
                # column.export_parameter_table.visible = True
                if column.operating_point_method.value == 'Grid':
                    column.nr_points_sampling.visible = True
                    column.consolidation_method.visible = True
                else:
                    column.nr_points_sampling.visible = False
                    column.consolidation_method.visible = False

    def update_column_block_neighbors(self, name, value):

        self.phenotypic_neighbors.user_cases.visible = False
        self.phenotypic_neighbors.all_mutation_rates_vbox.visible = False
        self.phenotypic_neighbors.phenotype.visible = True
        self.phenotypic_neighbors.export_table.visible = True
        self.phenotypic_neighbors.add_column.visible = True
        self.phenotypic_neighbors.add_filter.visible = True
        self.phenotypic_neighbors.action_button.description = 'Calculate'

        if self.phenotypic_neighbors.neighbors_method.value == 'Signature':
            self.phenotypic_neighbors.explicitly_share.visible = False
            self.phenotypic_neighbors.intersecting.visible = False
        elif self.phenotypic_neighbors.neighbors_method.value == 'Vertex Enumeration':
            self.phenotypic_neighbors.explicitly_share.visible = True
            self.phenotypic_neighbors.intersecting.visible = True
        elif self.phenotypic_neighbors.neighbors_method.value == 'User Defined':
            self.phenotypic_neighbors.explicitly_share.visible = False
            self.phenotypic_neighbors.intersecting.visible = False
            self.phenotypic_neighbors.user_cases.visible = True
        else:
            self.phenotypic_neighbors.phenotype.visible = False
            self.phenotypic_neighbors.all_mutation_rates_vbox.visible = True
            self.phenotypic_neighbors.export_table.visible = False
            self.phenotypic_neighbors.add_column.visible = False
            self.phenotypic_neighbors.add_filter.visible = False
            self.phenotypic_neighbors.intersecting.visible = False
            self.phenotypic_neighbors.explicitly_share.visible = False
            self.phenotypic_neighbors.action_button.description = 'Get Mutation Rates'

    def update_all_mutations_block(self, name, value):

        block = self.phenotypic_neighbors.all_mutation_rates

        block['asymmetrical_bounds_box'].visible = block['asymmetrical_bounds_checkbox'].value
        block['lower_bound_par_space'].visible = not (block['asymmetrical_bounds_checkbox'].value)
        block['upper_bound_par_space'].visible = not (block['asymmetrical_bounds_checkbox'].value)

        if block['nr_mutations'].value == 'Single':

                    for i in ['volume_calculation_method',
                              'tracks_calculation_method',
                              'include_fractional_volume_data']:
                        block[i].visible = True

                    for i in ['enforce_bounds_operating_point',
                              'enforce_bounds_volume',
                              'operating_point_calculation_method']:
                        block[i].visible = False

                    block['method_for_width_determination'].visible = True if block['volume_calculation_method'].value == 'Track Method' else False  # 4 and 2
                    block['include_real_volumes'].visible = True if block['include_fractional_volume_data'].value is True else False  #12 and 11

                    volume_method_widget = block['volume_calculation_method']                           #2
                    methods = ['Vertex Enumeration', 'Tolerances', 'Bounding Box', 'Geometric Mean T. & BB.', 'Track Method']
                    methods_dic = dict((i, i) for i in methods)
                    volume_method_widget.values = methods_dic
        else:
                    for i in ['operating_point_calculation_method',
                              'volume_calculation_method',
                              'enforce_bounds_operating_point',
                              'enforce_bounds_volume']:
                        block[i].visible = True
                    for i in ['tracks_calculation_method',
                              'method_for_width_determination',
                              'include_fractional_volume_data',
                              'include_real_volumes',
                              ]:
                        block[i].visible = False

                    volume_method_widget = block['volume_calculation_method']           #2
                    methods = ['Vertex Enumeration', 'Tolerances', 'Bounding Box', 'Geometric Mean T. & BB.']
                    methods_dic = dict((i, i) for i in methods)
                    volume_method_widget.values = methods_dic


def average_width(case, range_slice, pvals, track_parameter, orthogonal_parameter, orthogonal_bounds):

    pvals_copy = pvals.copy()
    av_width = 0
    for orthogonal_bound in orthogonal_bounds:
        pvals_copy[orthogonal_parameter] = 10**orthogonal_bound
        tolerance = case.vertices_1D_slice(pvals_copy,
                                           track_parameter,
                                           range_slice=range_slice,
                                           log_out=True)
        if len(tolerance) != 0:
            av_width += (max(tolerance)[0] - min(tolerance)[0])/2
        else:
            av_width += 0
    return av_width


def create_string_for_table(c):

    s = '<tr align=center><td style="padding:0 15px 0 15px;">{0}</td>' \
         '<td style="padding:0 15px 0 15px;">{1}</td>'.format(
        c.case_number + '&#1007' if (c.is_cyclical is False and c.is_unstable is True) else c.case_number,
        c.signature)

    return s


def round_sig(x, sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)


def mutation_key_to_case(case1, case2, symbol=None):
    key = 'k'

    try:
        if int(case1) < 10:
            key += '0'+ str(case1)
        else:
            key += str(case1)
    except:
        key += str(case1)

    if symbol is not None:
        key += str(symbol)

    try:
        if int(case2) < 10:
            key += '0' + str(case2)
        else:
            key += str(case2)
    except:
        key += str(case2)

    return key


def get_track_dictionary(controller, valid_cases, lb, ub, p_bounds, method='Vertex Enumeration', constraints=''):

    l_b = controller.pvals.copy()
    u_b = controller.pvals.copy()
    for key in l_b.keys():
        l_b[key] = lb[key]
        u_b[key] = ub[key]
    track_dictionary = dict((key, []) for key in l_b.keys())

    if method == 'Vertex Enumeration':
        # Let's calculate vertices for all phenotypes.
        limit_vertices = False
        maxVertices = 0
        for case_nr in valid_cases:
            case = controller.ds(case_nr) if constraints == '' else controller.ds(case_nr, constraints=constraints)
            volume, nr_vertices, vertices_matrix, operating_point = case.volume_lrs(l_b,
                                                                                    u_b,
                                                                                    maxVertices,
                                                                                    limit_vertices,
                                                                                    return_vertices_matrix=True)
            for vertex in range(int(nr_vertices)):
                for nr, key in enumerate(l_b.keys()):
                    track_dictionary[key].append(round(vertices_matrix[vertex][nr], 4))

    elif method == 'Bounding Box':
        for case_nr in valid_cases:
            case = controller.ds(case_nr) if constraints == '' else controller.ds(case_nr, constraints=constraints)
            bb = case.bounding_box(p_bounds=p_bounds, log_out=True)
            for key in bb.keys():
                if len(bb[key]) < 2:
                    continue
                track_dictionary[key].append(round(bb[key][0], 11))
                track_dictionary[key].append(round(bb[key][1], 11))

    for key in track_dictionary.keys():
        track_dictionary[key] = np.unique(track_dictionary[key])
        track_dictionary[key].sort()

    return track_dictionary

def get_orthogonal_box(track_variable, orthogonal_track_dictionary, exclude_track_variable = True):

    # retunrs a list of dictionaries containing boundaries for each box:
    # e.g. [{'KP': (-3.0, -2.0), 'bp': (-3.0, -2.0)}, {'KP': (-3.0, -2.0), 'bp': (-2.0, 1.0)}, ...]

    variables = orthogonal_track_dictionary.keys()
    if exclude_track_variable is True:
        orthogonal_variables = [var for var in variables if var != track_variable]
    else:
        orthogonal_variables = variables
    all_bounds = []

    for var in orthogonal_variables:
        all_bounds.append(orthogonal_track_dictionary[var])

    comb = list(itertools.product(*all_bounds))

    # convert into lists of dictionaries. One dictionary for each box.
    all_bounds_dic_list = []
    for box in comb:
        box_dir = {}
        for n, variable in enumerate(orthogonal_variables):
            box_dir[variable] = box[n]
        all_bounds_dic_list.append(box_dir)

    return all_bounds_dic_list


def process_orthogonal_track_dictionary(orthogonal_track_dictionary):

    # transform into pairs defining bounds
    orthogonal_track_pair_dictionary = {}

    for key in orthogonal_track_dictionary.keys():
        bounds = []
        for n, lb in enumerate(orthogonal_track_dictionary[key]):
            if n + 1 < len(orthogonal_track_dictionary[key]):
                bounds.append((lb, orthogonal_track_dictionary[key][n+1]))

        orthogonal_track_pair_dictionary[key] = bounds

    return orthogonal_track_pair_dictionary


def get_constraint_list(orthogonal_box):

    constraints = []
    vol_box = 1

    for orthogonal_variable in orthogonal_box.keys():
        lb_ = "%.18f" % (10 ** orthogonal_box[orthogonal_variable][0])
        ub_ = "%.18f" % (10 ** orthogonal_box[orthogonal_variable][1])
        constraints.append(orthogonal_variable + ' > ' + str(lb_))
        constraints.append(orthogonal_variable + ' < ' + str(ub_))
        vol_box *= orthogonal_box[orthogonal_variable][1] - orthogonal_box[orthogonal_variable][0]

    return constraints, vol_box


def get_additional_boxes(track_variable, orthogonal_boxes, valid_cases, p_bounds, controller, lb, ub,
                         method='Bounding Box'):

    expanded_orthogonal_boxes = []
    for orthogonal_box in orthogonal_boxes:
        new_boxes = get_additional_boxes_for_box(orthogonal_box,
                                                 track_variable,
                                                 valid_cases,
                                                 p_bounds,
                                                 controller,
                                                 lb, ub, method=method)
        for box in new_boxes:
            expanded_orthogonal_boxes.append(box)

    if len(expanded_orthogonal_boxes) == len(orthogonal_boxes):
        return_boxes = expanded_orthogonal_boxes
    else:
        return_boxes = get_additional_boxes(track_variable, expanded_orthogonal_boxes, valid_cases, p_bounds,
                                            controller, lb, ub, method=method)

    return return_boxes


def get_additional_boxes_for_box(orthogonal_box, track_variable, valid_cases, p_bounds, controller, lb, ub,
                                 method='Bounding Box'):

    constraints, vol_box = get_constraint_list(orthogonal_box)
    valid_cases_track = []

    for case_str in valid_cases:
        case = controller.ds(case_str, constraints=constraints)
        if case.is_valid(p_bounds=p_bounds) is True:
            valid_cases_track.append(case_str)

    orthogonal_track_dictionary = get_track_dictionary(controller, valid_cases_track, lb, ub, p_bounds,
                                                       method=method, constraints=constraints)
    orthogonal_track_pair_dictionary = process_orthogonal_track_dictionary(orthogonal_track_dictionary)
    orthogonal_boxes = get_orthogonal_box(track_variable, orthogonal_track_pair_dictionary)

    return orthogonal_boxes


def get_boundaries_for_track_variable(controller, valid_cases_track, lb, ub, p_bounds, track_variable, constraints='',
                                      method='Vertex Enumeration'):

    # This function should return a list of indices indicating the additional boundaries. Use parts of the functions:
    # get_track_dictionary() and process_orthogonal_track_dictionary()

    track_dictionary = get_track_dictionary(controller, valid_cases_track, lb, ub, p_bounds,
                                            method=method, constraints=constraints)
    orthogonal_track_pair_dictionary = process_orthogonal_track_dictionary(track_dictionary)
    boundaries_for_track_variable = orthogonal_track_pair_dictionary[track_variable]

    return boundaries_for_track_variable


def get_track_volume_for_case(controller, case_nr, track_variable, constraints, range_slice,
                              method='Average of Methods', p_bounds=None):

    case = controller.ds(case_nr, constraints=constraints)
    bb = case.bounding_box(p_bounds=p_bounds, log_out=True)
    pvals = case.valid_interior_parameter_set(p_bounds=p_bounds)

    height = 1
    height_av = 0

    for orthogonal_variable in bb.keys():
        if orthogonal_variable == track_variable:
            continue
        if len(bb[orthogonal_variable]) == 0 or len(bb[orthogonal_variable]) == 1:
            return 0
        height_i = max(bb[orthogonal_variable]) - min(bb[orthogonal_variable])
        height_av += height_i
        height *= height_i
        orthogonal_coordinate = min(bb[orthogonal_variable]) + height_i * 0.5
        pvals[orthogonal_variable] = 10 ** orthogonal_coordinate
    height_av = height_av/(len(bb.keys()) - 1)

    if method == 'Average of Methods':

        average_width_, shape = get_width_track_method(case, pvals, track_variable, bb, range_slice,
                                                       method='Average Width')
        centroid_width_, shape = get_width_track_method(case, pvals, track_variable, bb, range_slice,
                                                        method='Width at Orthogonal Centroid')

        vol_average_width = average_width_ * height
        vol_centroid_width = centroid_width_ * height
        track_volume = (vol_average_width + vol_centroid_width) / 2

    else:

        width, shape = get_width_track_method(case, pvals, track_variable, bb, range_slice, method=method)
        if method == 'Custom Formulas' and shape == 'Case 2':
            short = width[0]
            lon = width[1]
            track_volume = height * (0.5 * lon * height_av ** 2 - (lon - short) * (1.0/6.0) * height_av ** 3)
        elif method == 'Custom Formulas' and shape == 'Case 3':
            short = width[0]
            lon = width[1]
            track_volume = height * (0.5 * short * height_av ** 2 + (lon - short) * (1.0/6.0) * height_av ** 3)
        elif method == 'Custom Formulas' and shape == 'Case 4':
            track_volume = height * width / 3
        elif method == 'Custom Formulas' and shape == 'Case 5':
            track_volume = 0.5 * height * width * (1 + (height_av ** 2) / 2.0 + (height_av ** 3) / 6.0)
        elif method == 'Custom Formulas' and shape == 'Case 6':
            short = width[0]
            lon = width[1]
            slope = width[2]
            l_s = lon - short
            track_volume = ((2.0*lon)/2.0 - l_s/slope * (0.5 * l_s ** 2 - 1/3.0 * l_s ** 3 ) ) * height if slope != 0 else 0

        else:               # Case 1
            track_volume = height * width

    return track_volume, shape


def get_width_track_method(case, pvals_orthogonal, track_variable, bb, range_slice, method='Average Width'):

    shape = 'Case 1'
    width = 0
    pvals_copy = pvals_orthogonal.copy()

    if method == 'Average Width':

            orthogonal_vars = []
            all_bounds = []
            for key in bb.keys():
                if key == track_variable:
                    continue
                orthogonal_vars.append(key)
                all_bounds.append(bb[key])
            combs = list(itertools.product(*all_bounds))

            width = 0
            count = 0
            for comb in combs:

                for n, var in enumerate(orthogonal_vars):
                    pvals_copy[var] = 10 ** comb[n]

                tol = case.vertices_1D_slice(pvals_copy,
                                             track_variable,
                                             range_slice=range_slice[track_variable],
                                             log_out=True)
                if len(tol) != 0:
                    width += max(tol)[0] - min(tol)[0]
                    count += 1

            width = width/count if count != 0 else 0

    elif method == 'Width at Orthogonal Centroid':

            tol = case.vertices_1D_slice(pvals_orthogonal,
                                         track_variable,
                                         range_slice=range_slice[track_variable],
                                         log_out=True)
            if len(tol) == 0:
                width = 0
            else:
                width = max(tol)[0] - min(tol)[0]
            pvals_orthogonal[track_variable] = 10 ** (min(tol)[0] + 0.5 * width)

    elif method == 'Average Width over Facets at Centroid':

        # p_vals_orthogonal is already at the centroid.
        width = 0
        for orthogonal_var in bb.keys():
            if orthogonal_var == track_variable:
                continue

            width_i = 0
            pvals_copy = pvals_orthogonal.copy()
            for vertex in bb[orthogonal_var]:
                pvals_copy[orthogonal_var] = 10**vertex
                tol = case.vertices_1D_slice(pvals_copy,
                                             track_variable,
                                             range_slice=range_slice[track_variable],
                                             log_out=True)
                if len(tol) != 0:
                    width_i += max(tol)[0] - min(tol)[0]

            width += width_i/2
        width = width/(len(bb.keys()) - 1)

    elif method == 'Custom Formulas':

        orthogonal_vars = []
        all_bounds = []
        for key in bb.keys():
            if key == track_variable:
                continue
            orthogonal_vars.append(key)
            all_bounds.append(bb[key])
        combs = list(itertools.product(*all_bounds))

        width = 0
        count = 0
        full = 0
        dot = 0
        width_vector = []
        coordinates_vector = []
        coordinates_vector_max = []
        coordinates_vector_min = []
        missing = 0
        facet_number_vector = []
        n_facet = 1

        for jj, comb in enumerate(combs):

            for n, var in enumerate(orthogonal_vars):
                pvals_copy[var] = 10 ** comb[n]

            tol = case.vertices_1D_slice(pvals_copy,
                                         track_variable,
                                         range_slice=range_slice[track_variable],
                                         log_out=True)

            if len(tol) != 0:
                width += max(tol)[0] - min(tol)[0]
                width_vector.append(round(max(tol)[0] - min(tol)[0], 4))
                coordinates_vector.append([min(tol)[0], max(tol)[0]])
                coordinates_vector_max.append(max(tol)[0])
                coordinates_vector_min.append(min(tol)[0])
                facet_number_vector.append(n_facet)
                count += 1
            else:
                missing += 1

            if len(tol) == 1:
                dot += 1
            if len(tol) == 2:
                full += 1

            # if (jj + 1) % len(orthogonal_vars) == 0:
            n_facet += 1

        identical, different, zero = analyze_width_vector(width_vector)

        if full >= 2:

            if dot == 1 and identical >= 2:

                width = width/(count - dot)    # builds a non-zero average of lengths
                shape = 'Case 5'    #Spear Shortened
                return width, shape

            if full == 4 and identical == 3:

                shape = 'Case 6'
                facet_number_vector = np.array(facet_number_vector)
                coordinates_vector = np.array(coordinates_vector)
                width_min = min(width_vector)
                width_max = max(width_vector)
                coordinates_vector_max = np.array(coordinates_vector_max)
                mask = coordinates_vector_max == max(coordinates_vector_max)
                f = facet_number_vector[mask]
                if len(f) > 1:
                    f = f[0]
                delta = 2 if f <= 2 else - 2
                # print("facet_number_vector = ", facet_number_vector)
                # print("f, f+delta : ", f, f + delta)
                a_coordinates = coordinates_vector[facet_number_vector == f][0]
                a = min(a_coordinates)
                b_coordinates = coordinates_vector[facet_number_vector == f + delta][0]
                b = min(b_coordinates)
                slope = max(a, b) - min(a, b)

                if track_variable == 'KN' and case.case_number == '7':
                    print("short ", width_min)
                    print("long ", width_max)
                    print("slope ", slope)

                return [width_min, width_max, slope], shape

            if dot == 0 and identical == 2 and missing >= 1:
                unique, counts = np.unique(np.array(width_vector), return_counts=True)
                duplicates = unique[counts > 1]
                if max(width_vector) in duplicates:
                    shape = 'Case 2'  #'Spear Substracted'
                    width_min = min(width_vector)
                    width_max = max(width_vector)
                if min(width_vector) in duplicates:
                    shape = 'Case 3'
                    width_min = min(width_vector)
                    width_max = max(width_vector)
                return [width_min, width_max], shape

            width, shape = get_width_track_method(case, pvals_orthogonal, track_variable, bb, range_slice,
                                                  method='Average Width over Facets at Centroid')
            return width, shape
        else:
            shape = 'Case 4' # Pyramid
            width = width/2  ## update this

    return width, shape


def get_track_volume_for_case_decomposing(controller, case_nr, track_variable, constraints, range_slice,
                                          boundaries_for_track_variable, method='Average of Methods',
                                          p_bounds=None,
                                          decomposing=True,
                                          print_bol=False,
                                          tracks_method='Bounding Box'):
    total_vol = 0
    shape = 'None'
    if decomposing is False:
        total_vol, shape = get_track_volume_for_case(controller, case_nr, track_variable, constraints, range_slice,
                                                    method=method, p_bounds=p_bounds)
        return total_vol, shape


    boxes_for_polytope = get_boxes_for_polytope(controller, case_nr, track_variable, constraints, range_slice,
                                                boundaries_for_track_variable,
                                                tracks_method='Bounding Box', p_bounds=p_bounds)
    for band in boxes_for_polytope:
        constraints, vol_box = get_constraint_list(band)
        vol, shape = get_track_volume_for_case(controller, case_nr, track_variable, constraints, range_slice,
                                               method=method, p_bounds=p_bounds)
        total_vol += vol

    return total_vol, shape


def get_boxes_for_polytope(controller, case_nr, track_variable, constraints, range_slice,
                           boundaries_for_track_variable,
                           tracks_method='Bounding Box',
                           p_bounds=None):

    l_b = controller.pvals.copy()
    u_b = controller.pvals.copy()
    for key in l_b.keys():
        l_b[key] = min(range_slice)
        u_b[key] = max(range_slice)

    track_dictionary = dict((key, np.array([])) for key in l_b.keys())

    for band in boundaries_for_track_variable:

        additional_constraints, vol_box = get_constraint_list({str(track_variable): band})
        new_const = additional_constraints + constraints
        case_i = controller.ds(case_nr, constraints=new_const)

        if tracks_method == 'Vertex Enumeration':
            limit_vertices = False
            maxVertices = 0
            try:
                volume, nr_vertices, vertices_matrix, operating_point = case_i.volume_lrs(l_b,
                                                                                          u_b,
                                                                                          maxVertices,
                                                                                          limit_vertices,
                                                                                          return_vertices_matrix=True)
            except:
                continue
            for vertex in range(int(nr_vertices)):
                for nr, key in enumerate(l_b.keys()):
                    track_dictionary[key] = np.append(track_dictionary[key], round(vertices_matrix[vertex][nr], 4))

        elif tracks_method == 'Bounding Box':
            bb = case_i.bounding_box(p_bounds=p_bounds, log_out=True)

            for key in bb.keys():
                if len(bb[key]) < 2:
                    continue
                track_dictionary[key] = np.append(track_dictionary[key], round(bb[key][0], 11))
                track_dictionary[key] = np.append(track_dictionary[key], round(bb[key][1], 11))

        for key in track_dictionary.keys():
            track_dictionary[key] = np.unique(track_dictionary[key])
            track_dictionary[key].sort()

    ## Now process the track dictionary
    orthogonal_track_pair_dic = process_orthogonal_track_dictionary(track_dictionary)
    orthogonal_boxes = get_orthogonal_box(track_variable, orthogonal_track_pair_dic, exclude_track_variable=False)

    # orthogonal_boxes_out = orthogonal_boxes
    orthogonal_boxes_out = get_track_dictionary_for_polytope_recurrent(controller, case_nr, track_variable, constraints,
                                                                       l_b, u_b, orthogonal_boxes,
                                                                       tracks_method=tracks_method, p_bounds=p_bounds)

    return orthogonal_boxes_out


def get_track_dictionary_for_polytope_recurrent(controller, case_nr, track_variable, constraints, lb, ub,
                                                orthogonal_boxes_in,
                                                tracks_method='Bounding Box',
                                                p_bounds=None):

    track_dictionary = dict((key, np.array([])) for key in lb.keys())

    for box in orthogonal_boxes_in:

        additional_constraints, vol_box = get_constraint_list(box)
        new_const = additional_constraints + constraints
        case_i = controller.ds(case_nr, constraints=new_const)
        if tracks_method == 'Vertex Enumeration':
            limit_vertices = False
            maxVertices = 0
            try:
                volume, nr_vertices, vertices_matrix, operating_point = case_i.volume_lrs(lb,
                                                                                          ub,
                                                                                          maxVertices,
                                                                                          limit_vertices,
                                                                                          return_vertices_matrix=True)
            except:
                continue
            for vertex in range(int(nr_vertices)):
                for nr, key in enumerate(lb.keys()):
                    track_dictionary[key] = np.append(track_dictionary[key], round(vertices_matrix[vertex][nr], 4))

        elif tracks_method == 'Bounding Box':
            bb = case_i.bounding_box(p_bounds=p_bounds, log_out=True)

            for key in bb.keys():
                if len(bb[key]) < 2:
                    continue
                track_dictionary[key] = np.append(track_dictionary[key], round(bb[key][0], 11))
                track_dictionary[key] = np.append(track_dictionary[key], round(bb[key][1], 11))

        for key in track_dictionary.keys():
            track_dictionary[key] = np.unique(track_dictionary[key])
            track_dictionary[key].sort()

    ## Now process the track dictionary
    orthogonal_track_pair_dic = process_orthogonal_track_dictionary(track_dictionary)
    orthogonal_boxes_out = get_orthogonal_box(track_variable, orthogonal_track_pair_dic, exclude_track_variable=False)

    if len(orthogonal_boxes_in) == len(orthogonal_boxes_out) or len(orthogonal_boxes_out) > 100:
        return_boxes = orthogonal_boxes_out

    else:
        return_boxes = get_track_dictionary_for_polytope_recurrent(controller, case_nr, track_variable, constraints, lb,
                                                                   ub,
                                                                   orthogonal_boxes_out,
                                                                   tracks_method=tracks_method,
                                                                   p_bounds=None)
    return return_boxes


def is_height_zero(orthogonal_box, bounding_box_dic, from_case_key, to_case_key, pvals_from_case, pvals_to_case):

    height_zero = False
    height = 1

    # We calculate the height for the orthogonal parameters
    for orthogonal_variable in orthogonal_box.keys():
        height_i = min(bounding_box_dic[from_case_key][orthogonal_variable][1],
                       bounding_box_dic[to_case_key][orthogonal_variable][1]) - \
                   max(bounding_box_dic[from_case_key][orthogonal_variable][0],
                       bounding_box_dic[to_case_key][orthogonal_variable][0])

        # If height is negative or zero, then go to the next track parameter
        if height_i <= 0:
            height_zero = True
            continue

        height = height * height_i

        orthogonal_bounds = [max(bounding_box_dic[from_case_key][orthogonal_variable][0],
                                 bounding_box_dic[to_case_key][orthogonal_variable][0]),
                             min(bounding_box_dic[from_case_key][orthogonal_variable][1],
                                 bounding_box_dic[to_case_key][orthogonal_variable][1]),
                             ]

        orthogonal_coordinate = orthogonal_bounds[0] + height_i * 0.5
        pvals_from_case[orthogonal_variable] = 10 ** orthogonal_coordinate
        pvals_to_case[orthogonal_variable] = 10 ** orthogonal_coordinate

    return height_zero, height


def analyze_width_vector(width_vector):

    width_vector = np.array(width_vector)
    zero = sum(width_vector == 0)

    unique, counts = np.unique(width_vector, return_counts=True)
    different = len(unique)
    identical = sum(counts[counts > 1])



    # zero = 0
    # for element in width_vector:
    #     if element == 0:
    #         zero += 1

    # a = np.array(width_vector)
    # b = np.unique(a)
    # identical = len(a) - len(b) + 1
    # different = len(b) - 1

    return identical, different, zero


