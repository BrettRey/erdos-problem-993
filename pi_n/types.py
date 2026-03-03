from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional, Tuple, Union

UID = str


@dataclass(frozen=True, order=True)
class ABClass:
    a: int
    b: int

    def __post_init__(self) -> None:
        if self.a < 0 or self.b < 0:
            raise ValueError("a,b must be nonnegative")
        if self.a > self.b:
            raise ValueError("require a <= b")


@dataclass(frozen=True)
class InstanceRecord:
    uid: UID
    cls: ABClass
    Lambda: int
    D: int
    R_shift: int
    sum_all: int

    def __post_init__(self) -> None:
        if self.R_shift < 0:
            raise ValueError("R_shift must be nonnegative")


@dataclass(frozen=True)
class Ratio:
    num: int
    den: int

    def __post_init__(self) -> None:
        if self.den <= 0:
            raise ValueError("Ratio.den must be > 0")


Bucket = Tuple[UID, ...]
Partition = Tuple[Bucket, ...]


@dataclass(frozen=True)
class BucketCert:
    alpha: Ratio
    lambda_: Ratio
    witness_lower: UID
    witness_upper: UID


@dataclass(frozen=True)
class SplitStep:
    bucket_index: int
    cls: ABClass


@dataclass(frozen=True)
class Transcript:
    base: Literal["LINE_SPLIT_A_LE_3"]
    splits: Tuple[SplitStep, ...]


@dataclass(frozen=True)
class PiPass:
    status: Literal["PASS"]
    partition: Partition
    bucket_certs: Tuple[BucketCert, ...]
    transcript: Transcript
    input_digest: Optional[str] = None


@dataclass(frozen=True)
class DegenerateFailure:
    uid: UID
    cls: ABClass
    Lambda: int
    D: int
    sum_all: int
    violations: Tuple[Literal["sum_all > D", "Lambda < D - sum_all"], ...]


@dataclass(frozen=True)
class PiFailDegenerate:
    status: Literal["FAIL_DEGENERATE"]
    failures: Tuple[DegenerateFailure, ...]
    input_digest: Optional[str] = None


@dataclass(frozen=True)
class PiFailClass:
    status: Literal["FAIL_CLASS"]
    obstruction_class: ABClass
    obstruction_bucket_index: int
    bucket_cert: BucketCert
    transcript: Transcript
    input_digest: Optional[str] = None


PiResult = Union[PiPass, PiFailDegenerate, PiFailClass]
