// Exact row lower HNFs. Columns of H^T are a sublattice basis.
import { Q, matVec } from "./rational";
import { type IMat } from "./intmat";
export function indexKSublattices(n: number, k: number): IMat[] {
    if (!Number.isSafeInteger(n) || !Number.isSafeInteger(k) || n < 1 || k < 1)
        throw new RangeError("positive dimension and index required");
    const out: IMat[] = [];
    function diagonals(prefix: number[], left: number): void {
        if (prefix.length === n - 1) {
            fill([...prefix, left]);
            return;
        }
        for (let a = 1; a <= left; a++)
            if (left % a === 0)
                diagonals([...prefix, a], left / a);
    }
    function fill(diag: number[]): void {
        const H = diag.map((x, i) => diag.map((_, j) => BigInt(i === j ? x : 0)));
        const places = diag.flatMap((_, i) => Array.from({ length: i }, (_, j) => [i, j]));
        function visit(p: number): void {
            if (p === places.length) {
                out.push(H.map(r => [...r]));
                return;
            }
            const [i, j] = places[p];
            for (let v = 0; v < diag[j]; v++) {
                H[i][j] = BigInt(v);
                visit(p + 1);
            }
        }
        visit(0);
    }
    diagonals([], k);
    return out;
}
export function reduce(v: Q[], H: IMat): Q[] {
    let w = [...v];
    for (let i = H.length - 1; i >= 0; i--) {
        const q = w[i].div(Q.of(H[i][i])).floor();
        w = w.map((x, j) => x.sub(Q.of(q * H[i][j])));
    }
    return w;
}
export function isInvariant(H: IMat, Ws: IMat[]): boolean {
    return Ws.every(W => H.every(h => reduce(matVec(W.map(r => r.map(x => Q.of(x))), h.map(x => Q.of(x))), H).every(x => x.isZero())));
}
