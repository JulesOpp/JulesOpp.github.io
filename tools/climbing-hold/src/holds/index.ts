import type { HoldDef } from "../spec";
import { PARAMETRIC } from "./parametric";
import { SCULPTED } from "./sculpted";
import { MACHINABLE_EDGE } from "./machinable";

/** Every built-in hold, in menu order. Scanned holds are added at runtime. */
export const HOLDS: HoldDef[] = [MACHINABLE_EDGE, ...PARAMETRIC, ...SCULPTED];
export const holdById = (id: string) => HOLDS.find((h) => h.id === id);
