class DinkumException(Exception):
    pass

class DinkumObservationFailed(DinkumException):
    pass

class DinkumInvalidGene(DinkumException):
    pass

class DinkumInvalidTissue(DinkumException):
    pass

class DinkumNoSuchGene(DinkumException):
    pass

class DinkumInvalidActivationFunction(DinkumException):
    pass

class DinkumInvalidActivationResult(DinkumException):
    pass
