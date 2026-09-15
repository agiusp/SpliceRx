"""Data-load integration layer.

Reads a TCGA cohort directory that already sits on the server's disk (produced
offline by the ``prepTCGAdata`` R package), classifies the files it finds, and
loads the ones the apps can use straight into a session: the SJV junction matrix
+ sample metadata, the 2D View (``sjvc``) gene matrix + clinical table, and the
SJSurv sjdat matrices + sample metadata. It deliberately depends on the ``sjv``,
``sjvc`` and ``sjsurv`` packages (it is the glue between a cohort directory and
those apps).
"""
