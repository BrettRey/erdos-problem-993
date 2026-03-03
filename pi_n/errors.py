class PiContractError(Exception):
    """Base error for pi_n contract violations."""


class InputLogError(PiContractError):
    """Raised when the input log violates required invariants."""


class CertificateError(PiContractError):
    """Raised when a certificate fails replay verification."""
