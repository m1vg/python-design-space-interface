"""Box class.
Represents a container that can be used to group other widgets.
"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from ipywidgets import Widget, DOMWidget, register, widget_serialization, FlexBox
from traitlets import Unicode, Tuple, Int, CaselessStrEnum, Instance


class Popup(DOMWidget):
    """Displays multiple widgets in a group."""
    _model_name = Unicode('BoxModel', sync=True)
    _view_module = Unicode('nbextensions/popup/popup', sync=True)
    _view_name = Unicode('PopupView', sync=True)

    # Child widgets in the container.
    # Using a tuple here to force reassignment to update the list.
    # When a proper notifying-list trait exists, that is what should be used here.
    children = Tuple(sync=True, **widget_serialization)

    _overflow_values = ['visible', 'hidden', 'scroll', 'auto', 'initial', 'inherit', '']
    overflow_x = CaselessStrEnum(
        values=_overflow_values,
        default_value='auto', sync=True, help="""Specifies what
        happens to content that is too large for the rendered region.""")
    overflow_y = CaselessStrEnum(
        values=_overflow_values,
        default_value='auto', sync=True, help="""Specifies what
        happens to content that is too large for the rendered region.""")

    box_style = CaselessStrEnum(
        values=['success', 'info', 'warning', 'danger', ''],
        default_value='', allow_none=True, sync=True, help="""Use a
        predefined styling for the box.""")

    def __init__(self, children = (), **kwargs):
        kwargs['children'] = children
        super(Popup, self).__init__(**kwargs)
        self.on_displayed(Popup._fire_children_displayed)
        self.width='45%'
        self.height='45%'

    def _fire_children_displayed(self):
        for child in self.children:
            child._handle_displayed()

            