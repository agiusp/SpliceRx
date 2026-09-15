"""Standalone SJV FastAPI application.

Kept for running / testing SJV in isolation. The combined service is
``app/backend/main.py``, which mounts this package's router under ``/api/sjv``.
Here the router is served under a plain ``/api`` prefix.
"""
from __future__ import annotations

from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .api.routes import router
from .services.sessions import store


@asynccontextmanager
async def lifespan(app: FastAPI):
    store.start_sweeper()
    yield


app = FastAPI(title="SJV — Splice Junction Visualizer", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(router, prefix="/api")


@app.get("/api/health")
def health() -> dict:
    return {"ok": True}
