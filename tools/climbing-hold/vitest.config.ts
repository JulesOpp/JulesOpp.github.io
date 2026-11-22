import { defineConfig } from "vitest/config";

export default defineConfig({
  test: {
    include: ["test/**/*.test.ts"],
    // OpenCascade and manifold wasm are heavy; never use more than 4 workers
    maxWorkers: 4,
  },
});
