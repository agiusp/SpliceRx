import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// One dev server for both tabs. Everything under /api is proxied to the single
// combined backend, which routes /api/sjv/* and /api/sjvc/* internally.
export default defineConfig({
  plugins: [react()],
  server: {
    port: 5173,
    proxy: {
      "/api": "http://localhost:8000",
    },
  },
});
