"""
If you have the enthought python distribution, or install traits+traitsui
you can then examine your NMR glue dictionary, by writing

from nmrglue.util import dicgui
dicgui.dic(dic)
"""

from traits.api import HasTraits, Instance
from traitsui.api import View, VGroup, Item, ValueEditor

class DictEditor(HasTraits):
    Object = Instance( object )
    def __init__(self, obj, **traits):
        super(DictEditor, self).__init__(**traits)
        self.Object = obj
    def trait_view(self, name=None, view_elements=None):
        return View(
          VGroup(
            Item('Object',
                  label      = 'Debug',
                  id         = 'debug',
                  editor     =ValueEditor(), #ValueEditor()
                  style      = 'custom',
                  dock       = 'horizontal',
                  show_label = False),),
          title     = 'Dictionary Editor',
          width     = 800,
          height    = 600,
          resizable = True)
def dic(my_data):
    b = DictEditor(my_data)
    b.configure_traits()
