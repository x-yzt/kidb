import operator
from collections import OrderedDict
from collections.abc import MutableSequence
from functools import reduce

import pandas


def is_outlier(deviation=3):
    """Build a callable taking a `Series` as argument, and returning a new
    `Series` of booleans indicating whether values are in a specific SD range
    or not.
    """
    def func(series):
        mean, std = series.mean(), series.std()

        return ~series.between(
            mean - deviation * std,
            mean + deviation * std
        )
    return func


class KiDatabase:
    """Wrapper around a `DataFrame`, allowing it to be used somewhat like a
    common web database - as weird as it may sound.
    """
    structure = OrderedDict({
        "Name":         'receptor',
        "Unigene":      'unigene',
        " Ligand Name": 'ligand',
        "CAS":          'cas',
        "NSC":          'nsc',
        "Hotligand":    'ref_ligand',
        "species":      'species',
        "source":       'source',
        "ki Note":      'ki_op',
        "ki Val":       'ki',
        "Reference":    'reference',
        "Link":         'link',
    })

    def __init__(self, path):
        """Create a new class instance by loading values from a CSV file in a
        new DataFrame.
        """
        df = pandas.read_csv(path, usecols=self.structure.keys())
        df.rename(columns=self.structure, inplace=True)
        df['ki_op'].fillna(value='=', inplace=True)
        for col in df.columns:
            if col != 'ki':
                df[col].fillna('', inplace=True)

        self.data = df

    def __str__(self):
        return str(self.data)

    @classmethod
    def from_dataframe(cls, df):
        """Create a new class instance wrapping a bare DataFrame."""
        self = cls.__new__(cls)
        self.data = df
        return self

    @property
    def columns(self):
        return self.data.columns

    def filter(self, **kwargs):
        """Simple value filter. Return a new class instance.

        Multiple arguments will be combinated using the AND operator.

        Arguments have to be `key=value` pairs. Multiple values, such as
        `key=[value1, value2]`, are supported and will be combinated using the
        OR operator.
        """
        if not kwargs:
            return self

        for key, value in kwargs.items():
            if not isinstance(value, MutableSequence):
                kwargs[key] = [value]

        return type(self).from_dataframe(
            self.data.loc[
                reduce(
                    operator.and_,
                    (self.data[k].isin(v) for k, v in kwargs.items())
                )
            ]
        )

    def values(self, column):
        """Get all unique values from a given column."""
        return self.data[column].unique()

    def get_ki(self, ligand, deviation=None):
        """Get various aggregations of Ki data for a given ligand. All the
        sources used for calculations will also be returned.

        `deviation` is in standard deviation units.
        Note excluded outliers will not affect computed values.
        """
        df = self.filter(ligand=ligand, ki_op='=').data.sort_values('receptor')
        ki = df.groupby('receptor')['ki']

        if deviation and not df.empty:
            df = df[~ki.apply(is_outlier(deviation))]
            ki = df.groupby('receptor')['ki']

        return {
            'ki': pandas.DataFrame({
                'median': ki.median(),
                'mean':   ki.mean(),
                'std':    ki.std(),
                'count':  ki.count(),
            }),
            'sources': df
        }
