"""SJVC FastAPI application."""
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


app = FastAPI(title="SJVC — Splice Junction Visualizer with Clinical", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5174", "http://127.0.0.1:5174"],
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(router)
