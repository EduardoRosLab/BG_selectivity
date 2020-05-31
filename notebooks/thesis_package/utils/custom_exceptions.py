
class ParamError(Exception):
    """Any kind of error related to parameters values."""
    def __init__(self, param_name, msg=None):
        msg = "Some problem occured with parameter {}.".format(param_name) if msg is None else msg
        super().__init__(msg)
        self.param_name = param_name