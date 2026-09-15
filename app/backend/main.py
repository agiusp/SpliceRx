"""Combined SJ Web Apps service.

Mounts the previously separate apps behind one FastAPI process. Their routers
keep their own code and session stores; they are only separated by URL prefix:

    /api/sjv/*      SJV — sashimi visualizer
    /api/sjvc/*     SJVC — cohort heatmap / projection (a.k.a. 2D View in the UI)
    /api/sjsurv/*   SJSurv — survivor-group classification from splice junctions
    /api/sjlookup/* SJ Lookup — direct per-junction lookup against junction metadata
    /api/dataload/* Data tab — load a TCGA cohort directory into a session
    /api/health     combined health check

Run from this directory:

    uvicorn main:app --port 8000 --reload
"""
from __future__ import annotations

from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from dataload.routes import router as dataload_router
from sjlookup.api.routes import router as sjlookup_router
from sjlookup.services.sessions import store as sjlookup_store
from sjsurv.api.routes import router as sjsurv_router
from sjsurv.services.sessions import store as sjsurv_store
from sjv.api.routes import router as sjv_router
from sjv.services.sessions import store as sjv_store
from sjvc.api.routes import router as sjvc_router
from sjvc.services.sessions import store as sjvc_store

# The single Vite dev server that serves every tab.
FRONTEND_ORIGINS = ["http://localhost:5173", "http://127.0.0.1:5173"]


@asynccontextmanager
async def lifespan(app: FastAPI):
    sjv_store.start_sweeper()
    sjvc_store.start_sweeper()
    sjsurv_store.start_sweeper()
    sjlookup_store.start_sweeper()
    yield


app = FastAPI(title="SJ Web Apps", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=FRONTEND_ORIGINS,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(sjv_router, prefix="/api/sjv")
app.include_router(sjvc_router, prefix="/api/sjvc")
app.include_router(sjsurv_router, prefix="/api/sjsurv")
app.include_router(sjlookup_router, prefix="/api/sjlookup")
app.include_router(dataload_router, prefix="/api/dataload")


@app.get("/api/health")
def health() -> dict:
    return {"ok": True}
