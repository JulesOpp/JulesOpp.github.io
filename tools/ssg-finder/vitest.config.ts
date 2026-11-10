import { defineConfig } from "vitest/config";

export default defineConfig({
  test: {
    include: ["test/**/*.test.ts"],
    pool: "forks",
    testTimeout: 300_000,
    teardownTimeout: 5_000,
  },
});
