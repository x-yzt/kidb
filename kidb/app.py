from contextlib import suppress
from urllib.parse import unquote

from sanic import Sanic
from sanic.response import json

from kidb import BASE_DIR
from kidb.db import KiDatabase


app = Sanic("KiDB")
app.ctx.db = KiDatabase(BASE_DIR / "database.csv")


def get_from_params(request, param, type_=str, default=None):
    """Extract a value from a given request GET parameter."""
    value = None
    with suppress(ValueError):
        value = request.args.get(param)

    return default if value is None else type_(unquote(value))


# Middleware

@app.middleware('request')
async def extract_filters(request):
    request.ctx.filters = {
        k: v for k, v in request.args.items()
        if k in app.ctx.db.columns
    }


@app.middleware('response')
async def allow_cors(_, response):
    response.headers['Access-Control-Allow-Origin'] = '*'


# API routes

@app.route('/ligands', version=1)
async def ligands(request):
    max_len = get_from_params(request, 'max-len', int, -1)

    return json(list(
        filter(
            lambda x: len(x) <= max_len or max_len <= 0,
            app.ctx.db
                .filter(**request.ctx.filters)
                .values('ligand')
        )
    ))


@app.route('/receptors', version=1)
async def receptors(request):
    return json(
        app.ctx.db
            .filter(**request.ctx.filters)
            .values('receptor')
            .tolist()
    )


@app.route('/ki/<ligand:string>', version=1)
async def ki(request, ligand):
    ligand = unquote(ligand)
    deviation = get_from_params(request, 'deviation', float)

    res = (
        app.ctx.db
            .filter(**request.ctx.filters)
            .get_ki(ligand, deviation=deviation)
    )
    return json({
        k: (v
            .fillna('NaN')
            .to_dict(orient=('index', 'records')[k == 'sources'])
        )
        for k, v in res.items()
    })
