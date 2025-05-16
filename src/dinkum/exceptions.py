class DinkumException(Exception):
    pass


class DinkumObservationFailed(DinkumException):
    pass


class DinkumMultipleRules(DinkumException):
    pass


class DinkumInvalidGene(DinkumException):
    pass


class DinkumNotATranscriptionFactor(DinkumException):
    pass


class DinkumNotALigand(DinkumException):
    pass


class DinkumInvalidTissue(DinkumException):
    pass


class DinkumNoSuchGene(DinkumException):
    pass


class DinkumInvalidActivationFunction(DinkumException):
    pass


class DinkumInvalidActivationResult(DinkumException):
    pass
